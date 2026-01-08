import Pkg
# 自动检测并安装缺失的包
required_packages = ["DSP", "Optim"]
for pkg in required_packages
    if Base.find_package(pkg) === nothing
        println("正在安装 $pkg ...")
        Pkg.add(pkg)
    end
end

using DSP
using TyPlot
using Optim
using LinearAlgebra
using Statistics

# ==============================================================================
# 1. 辅助函数：构造全通滤波器
# ==============================================================================
# 根据极点半径(r)和角度(theta)构建级联的全通滤波器 (ZeroPoleGain形式)
function make_allpass(params)
    n_sections = div(length(params), 2)
    total_filter = nothing
    
    for i in 1:n_sections
        r = params[2*i - 1]
        theta = params[2*i]
        
        # 共轭极点对
        p = r * exp(im * theta)
        poles = [p, conj(p)]
        # 全通滤波器的零点是极点的倒数共轭
        zeros_vec = [1/conj(p), 1/p]
        
        # 增益设为 r^2 以保持单位增益
        section = ZeroPoleGain(zeros_vec, poles, 1.0)
        
        if total_filter === nothing
            total_filter = section
        else
            total_filter = total_filter * section 
        end
    end
    
    return total_filter
end

# ==============================================================================
# 2. 优化函数
# ==============================================================================
function optimize_group_delay(iir_filter, wp, N_ap)
    println("正在优化全通均衡器参数...")
    
    n_sections = div(N_ap, 2)
    # 优化频率范围：0 到 截止频率 (通带)
    w_pass = range(0, stop=wp*pi, length=100)
    
    # 计算原始 IIR 的群延时
    tau_iir = grpdelay(iir_filter, w_pass)
    mean_tau_iir = mean(tau_iir)
    
    # 目标函数：最小化 (总群延时 - 平均值) 的标准差
    function cost_function(params)
        try
            ap_filter = make_allpass(params)
            tau_ap = grpdelay(ap_filter, w_pass)
            tau_total = tau_iir + tau_ap
            # 我们希望群延时尽可能平坦，即标准差最小
            return std(tau_total)
        catch
            return Inf
        end
    end
    
    # 初始猜测 [r1, theta1, r2, theta2, ...]
    initial_params = Float64[]
    for i in 1:n_sections
        push!(initial_params, 0.5 + 0.1*i) # r
        push!(initial_params, 0.1*pi * i)  # theta
    end
    
    # 约束：半径 0~1 (稳定)，角度 0~pi
    lower_bounds = repeat([0.0, 0.0], n_sections)
    upper_bounds = repeat([0.99, pi], n_sections)
    
    # 使用 Fminbox + BFGS 优化
    res = optimize(cost_function, lower_bounds, upper_bounds, initial_params, Fminbox(BFGS()), Optim.Options(time_limit=15.0))
    
    best_params = Optim.minimizer(res)
    println("优化结束。最终代价 (std): ", Optim.minimum(res))
    return make_allpass(best_params)
end

# ==============================================================================
# 3. 绘图函数
# ==============================================================================
function plot_results(iir_filter, ap_filter, wp, ws)
    # --- 准备数据 ---
    # 频率轴：全频段 [0, pi] 用于幅频响应，通带 [0, wp*pi] 用于群延时细节
    n_points = 1000
    w_full = range(0, stop=pi, length=n_points)
    w_pass = range(0, stop=wp*pi, length=n_points)
    
    # 1. 计算幅频响应 (Magnitude)
    # 使用 freqz 计算复数响应
    H_iir = freqz(iir_filter, w_full)
    H_ap  = freqz(ap_filter, w_full)
    H_total = H_iir .* H_ap
    
    mag_iir = abs.(H_iir)
    mag_total = abs.(H_total)
    
    # 归一化 (以 DC 增益为基准)
    mag_iir ./= mag_iir[1]
    mag_total ./= mag_total[1] # 级联后可能会有微小增益变化，重新归一化
    
    # 2. 计算群延时 (Group Delay) - 仅关注通带和过渡带
    tau_iir = grpdelay(iir_filter, w_pass)
    tau_ap  = grpdelay(ap_filter, w_pass)
    tau_total = tau_iir + tau_ap
    
    # --- 开始绘图 ---
    figure("IIR Filter & Group Delay Equalization")
    clf()
    
    # 子图 1: 幅频响应 (Magnitude)
    subplot(2, 1, 1)
    hold("on")
    plot(w_full/pi, mag_iir, "b-", linewidth=1.5, label="原始低通 IIR")
    plot(w_full/pi, mag_total, "r--", linewidth=1.5, label="级联后 (IIR + 全通)")
    
    # 标记通带截止频率
    plot([wp, wp], [0, 1.1], "k:", linewidth=1)
    
    # 样式设置
    ylabel("归一化幅值")
    title("1. 幅频响应对比 (验证全通特性)")
    grid("on")
    legend()
    xlim([0, 1])
    ylim([0, 1.1])
    # 可以在图上标注一下，全通不改变幅度
    # 修复：将 color="gray" 改为十六进制 "#808080"，因为 TyPlot 不支持 "gray" 名称
    text(0.5, 0.5, "注意：红虚线应与蓝实线重合\n(全通滤波器不改变幅度谱)", fontsize=8, color="#808080")
    
    # 子图 2: 群延时 (Group Delay)
    subplot(2, 1, 2)
    hold("on")
    plot(w_pass/pi, tau_iir, "b-", linewidth=1.5, label="原始群延时")
    plot(w_pass/pi, tau_total, "r-", linewidth=2.0, label="均衡后总群延时")
    
    # 绘制平均延时参考线
    avg_delay = mean(tau_total)
    plot([0, wp], [avg_delay, avg_delay], "g-.", linewidth=1, label="平均延时 ($(round(avg_delay, digits=1)))")
    
    # 样式设置
    xlabel("归一化频率 (×π rad/sample)")
    ylabel("群延时 (samples)")
    title("2. 通带群延时均衡效果")
    grid("on")
    legend()
    xlim([0, wp]) # 聚焦通带
    
    # 自动调整Y轴范围以展示细节
    y_min = minimum(tau_total) * 0.8
    y_max = maximum(tau_iir) * 1.1
    ylim([y_min, y_max])
    
    hold("off")
end

# ==============================================================================
# 主程序
# ==============================================================================
println("=== 开始设计 ===")

# 1. 滤波器指标
wp = 0.4    # 通带截止 (x pi)
ws = 0.5    # 阻带起始 (x pi)
Rp = 0.6    # dB
Rs = 32.0   # dB
N_iir = 3   # IIR 阶数

# 2. 设计 IIR 低通滤波器
println("1. 设计 3阶 椭圆低通滤波器...")
# 显式指定 design_method
response = Lowpass(wp)
method = DSP.Elliptic(N_iir, Rp, Rs)
iir_filter = digitalfilter(response, method)

# 3. 设计全通均衡器
println("2. 设计 6阶 全通均衡器 (优化中)...")
N_ap = 6
ap_filter = optimize_group_delay(iir_filter, wp, N_ap)

# 4. 绘图
println("3. 绘制分析图...")
plot_results(iir_filter, ap_filter, wp, ws)

println("=== 完成 ===")