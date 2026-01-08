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
        
        # 增益设为 r^2 以保持单位增益 (简单起见，后续绘图归一化处理)
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
    println("正在优化 10 阶全通均衡器参数 (这可能需要一点时间)...")
    
    n_sections = div(N_ap, 2)
    # 优化频率范围：0 到 截止频率 (通带)
    w_pass = range(0, stop=wp*pi, length=150) # 增加点数提高精度
    
    # 计算原始 IIR 的群延时
    tau_iir = grpdelay(iir_filter, w_pass)
    
    # 目标函数：最小化 (总群延时) 的标准差
    function cost_function(params)
        try
            ap_filter = make_allpass(params)
            tau_ap = grpdelay(ap_filter, w_pass)
            tau_total = tau_iir + tau_ap
            # 我们希望群延时尽可能平坦
            return std(tau_total)
        catch
            return Inf
        end
    end
    
    # 初始猜测 [r1, theta1, r2, theta2, ...]
    # 针对 10 阶 (5个节)，我们需要更细致的初始分布
    initial_params = Float64[]
    for i in 1:n_sections
        # 半径分布在 0.5 ~ 0.9 之间
        push!(initial_params, 0.5 + 0.08*i) 
        # 角度均匀分布在通带内，略微偏移
        push!(initial_params, (i / (n_sections+1)) * wp * pi)
    end
    
    # 约束：半径 0~0.99 (稳定)，角度 0~pi
    lower_bounds = repeat([0.0, 0.0], n_sections)
    upper_bounds = repeat([0.99, pi], n_sections)
    
    # 使用 Fminbox + BFGS 优化
    res = optimize(cost_function, lower_bounds, upper_bounds, initial_params, Fminbox(BFGS()), Optim.Options(time_limit=30.0, iterations=1000))
    
    best_params = Optim.minimizer(res)
    println("优化结束。最终标准差 (std): ", round(Optim.minimum(res), digits=4))
    return make_allpass(best_params)
end

# ==============================================================================
# 3. 绘图函数
# ==============================================================================
function plot_results(iir_filter, ap_filter, wp, ws)
    # --- 准备数据 ---
    n_points = 1000
    w_full = range(0, stop=pi, length=n_points)
    w_pass = range(0, stop=wp*pi, length=n_points)
    
    # 1. 计算幅频响应 (Magnitude)
    H_iir = freqz(iir_filter, w_full)
    H_ap  = freqz(ap_filter, w_full)
    H_total = H_iir .* H_ap
    
    mag_iir = abs.(H_iir)
    mag_ap = abs.(H_ap)
    mag_total = abs.(H_total)
    
    # 归一化 (以 DC 增益为基准)
    mag_iir ./= mag_iir[1]
    # 全通滤波器理论上幅度为1，但也归一化一下以防万一
    mag_ap ./= mag_ap[1]
    mag_total ./= mag_total[1]
    
    # 2. 计算群延时 (Group Delay)
    tau_iir = grpdelay(iir_filter, w_pass)
    tau_ap  = grpdelay(ap_filter, w_pass)
    tau_total = tau_iir + tau_ap
    
    # --- 开始绘图 ---
    figure("Q12: IIR Filter & Group Delay Equalization")
    clf()
    
    # 子图 1: 幅频响应
    subplot(2, 1, 1)
    hold("on")
    plot(w_full/pi, mag_iir, "b-", linewidth=1.5, label="原始低通 IIR")
    # 绘制全通滤波器的幅度 (应该是平直的)
    plot(w_full/pi, mag_ap, "g-.", linewidth=1.0, label="全通滤波器 (幅度)")
    plot(w_full/pi, mag_total, "r--", linewidth=1.5, label="级联后总响应")
    
    # 标记通带截止
    plot([wp, wp], [0, 1.2], "k:", linewidth=1)
    
    ylabel("归一化幅值")
    title("1. 幅频响应 (验证全通特性)")
    grid("on")
    legend("loc", "best")
    xlim([0, 1])
    ylim([0, 1.2])
    
    # 子图 2: 群延时
    subplot(2, 1, 2)
    hold("on")
    plot(w_pass/pi, tau_iir, "b-", linewidth=1.5, label="原始群延时")
    plot(w_pass/pi, tau_total, "r-", linewidth=2.0, label="均衡后总群延时")
    
    # 平均延时参考线
    avg_delay = mean(tau_total)
    plot([0, wp], [avg_delay, avg_delay], "g-.", linewidth=1, label="平均延时 ($(round(avg_delay, digits=1)))")
    
    xlabel("归一化频率 (×π rad/sample)")
    ylabel("群延时 (samples)")
    title("2. 通带群延时均衡效果 (10阶全通)")
    grid("on")
    legend("loc", "best")
    xlim([0, wp])
    
    # 自动缩放 Y 轴，留出一点余量
    y_min = minimum(tau_total) * 0.8
    y_max = maximum(tau_iir) * 1.1
    ylim([y_min, y_max])
    
    hold("off")
end

# ==============================================================================
# 主程序
# ==============================================================================
println("=== 第12题：设计开始 ===")

# 1. 滤波器指标
wp = 0.35   # 通带截止 0.35pi
ws = 0.5    # 阻带起始 (假设值，题目只给了最小阻带衰减，通常配合椭圆设计自适应，或者我们可以给一个合理的过渡带)
            # 注意：DSP.Elliptic 不需要显式 ws，它根据 N, Rp, Rs 设计
Rp = 0.8    # dB
Rs = 35.0   # dB
N_iir = 5   # IIR 阶数

# 2. 设计 IIR 低通滤波器
println("1. 设计 5阶 椭圆低通滤波器...")
response = Lowpass(wp)
method = DSP.Elliptic(N_iir, Rp, Rs)
iir_filter = digitalfilter(response, method)

# 3. 设计全通均衡器
println("2. 设计 10阶 全通均衡器...")
N_ap = 10
ap_filter = optimize_group_delay(iir_filter, wp, N_ap)

# 4. 绘图
println("3. 绘制分析图...")
plot_results(iir_filter, ap_filter, wp, ws)

println("=== 完成 ===")