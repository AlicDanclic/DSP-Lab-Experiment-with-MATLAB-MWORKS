import Pkg
# 自动检测并安装缺失的包 (DSP, Optim)
# TyPlot 通常在 MWorks 中预装，但为了保险起见也可以检查，不过这里主要解决 Optim
required_packages = ["DSP", "Optim"]
for pkg in required_packages
    if Base.find_package(pkg) === nothing
        println("正在安装 $pkg ...")
        Pkg.add(pkg)
    end
end

using DSP
using TyPlot # 替换 Plots 为 TyPlot
using Optim
using LinearAlgebra
using Statistics # 确保 mean 和 std 可用

# ==========================================
# 1. 设计 3 阶椭圆低通滤波器 (IIR)
# ==========================================
println("正在设计椭圆滤波器...")

# 性能指标
wp = 0.4             # 通带截止频率 (归一化频率, 1.0对应π)
ws = 0.5             # 阻带起始频率
Rp = 0.6             # 通带波纹 (dB)
Rs = 32.0            # 最小阻带衰减 (dB)
N_iir = 3            # 滤波器阶数

# 设计滤波器对象
# 显式使用 DSP.Elliptic 以避免与 TyMath.Elliptic 冲突
response_type = Lowpass(wp)
design_method = DSP.Elliptic(N_iir, Rp, Rs)
iir_filter = digitalfilter(response_type, design_method)

# ==========================================
# 2. 设计 6 阶全通均衡器 (优化过程)
# ==========================================
println("正在优化全通均衡器参数 (可能需要几秒钟)...")

# 均衡器阶数 (6阶 = 3个二阶节)
N_ap = 6
n_sections = div(N_ap, 2) 

# 定义频率网格 (只关注通带内的延时平坦度)
# 使用 0 到 0.4π 之间的频率点进行优化
w_pass = range(0, stop=wp*pi, length=100)

# 计算原始 IIR 滤波器的群延时
tau_iir = grpdelay(iir_filter, w_pass)

# 辅助函数：根据极点半径(r)和角度(theta)构建全通滤波器对象
function make_allpass(params)
    total_filter = nothing
    for i in 1:n_sections
        r = params[2*i - 1]
        theta = params[2*i]
        
        p = r * exp(im * theta)
        poles = [p, conj(p)]
        zeros = [1/conj(p), 1/p]
        
        section = ZeroPoleGain(zeros, poles, 1.0)
        
        if total_filter === nothing
            total_filter = section
        else
            total_filter = total_filter * section 
        end
    end
    return total_filter
end

# 目标函数：计算总群延时的标准差 (越小越平坦)
function cost_function(params)
    try
        ap_filter = make_allpass(params)
        tau_ap = grpdelay(ap_filter, w_pass)
        tau_total = tau_iir + tau_ap
        return std(tau_total) 
    catch
        return Inf 
    end
end

# 初始猜测参数 [r, theta, r, theta, ...]
initial_params = [0.5, 0.1*pi, 0.6, 0.2*pi, 0.7, 0.3*pi]

# 设置参数边界
lower_bounds = repeat([0.0, 0.0], n_sections)
upper_bounds = repeat([0.99, pi], n_sections)

# 执行优化
res = optimize(cost_function, lower_bounds, upper_bounds, initial_params, Fminbox(BFGS()), Optim.Options(time_limit=10.0))
best_params = Optim.minimizer(res)

println("优化完成。")
println("最佳全通极点参数 (r, theta): ", round.(best_params, digits=3))

# ==========================================
# 3. 结果计算与绘图 (TyPlot 版本)
# ==========================================

# 构建最终的全通滤波器
ap_filter_final = make_allpass(best_params)

# 在更宽的频率范围内计算响应以便绘图
w_plot = range(0, stop=wp*pi, length=500) 
tau_iir_plot = grpdelay(iir_filter, w_plot)
tau_ap_plot = grpdelay(ap_filter_final, w_plot)
tau_total_plot = tau_iir_plot + tau_ap_plot

# 准备绘图数据
x_data = w_plot ./ pi # 归一化频率

# 计算Y轴范围以便绘制垂直线
y_min = minimum([minimum(tau_iir_plot), minimum(tau_total_plot)])
y_max = maximum([maximum(tau_iir_plot), maximum(tau_total_plot)])

# 使用 TyPlot 绘图
# TyPlot 通常支持类似 plot(x, y, spec) 的语法
println("正在绘制结果...")

# 绘制原始低通滤波器群延时
# "b" 代表蓝色
plot(x_data, tau_iir_plot, "b", label="Original Lowpass")
hold("on") # 保持图像以叠加

# 绘制均衡后的群延时
# "r" 代表红色
plot(x_data, tau_total_plot, "r", label="Equalized (Total)")

# 绘制通带截止频率的垂直虚线
# 手动构造线段数据: ([wp, wp], [y_min, y_max])
# "k--" 代表黑色虚线
plot([wp, wp], [y_min, y_max], "k--", label="Cutoff Frequency")

title("Group Delay Equalization Result")
xlabel("Normalized Frequency (x pi rad/sample)")
ylabel("Group Delay (samples)")
legend() # 显示图例
grid("on") # 显示网格

println("平均群延时 (原始): ", round(mean(tau_iir_plot), digits=2))
println("平均群延时 (均衡后): ", round(mean(tau_total_plot), digits=2))
println("群延时波动 (原始 std): ", round(std(tau_iir_plot), digits=4))
println("群延时波动 (均衡后 std): ", round(std(tau_total_plot), digits=4))