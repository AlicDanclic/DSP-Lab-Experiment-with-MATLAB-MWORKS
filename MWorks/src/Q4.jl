import Pkg
# 检查并安装 DSP 包
required_packages = ["DSP", "TyPlot"]
for pkg in required_packages
    if Base.find_package(pkg) === nothing
        println("正在安装 $pkg ...")
        Pkg.add(pkg)
    end
end

using DSP
using TyPlot
using Printf

println("=== 切比雪夫 I 型 IIR 高通滤波器设计 ===")

# ==========================================
# 1. 定义性能指标
# ==========================================
Rp = 0.5            # 通带波纹 (dB)
Rs = 43.0           # 最小阻带衰减 (dB)
wp = 0.75           # 通带边缘频率 (归一化, 1.0 = π)
ws = 0.5           # 阻带边缘频率 (归一化, 1.0 = π)

println("指标参数:")
println("  通带波纹 Rp: $Rp dB")
println("  阻带衰减 Rs: $Rs dB")
println("  通带频率 wp: $(wp)π")
println("  阻带频率 ws: $(ws)π")

# ==========================================
# 2. 计算最小阶数 N
# ==========================================
# 为了确定满足衰减指标的最小阶数，我们需要使用模拟原型的公式
# 步骤: 
# 1. 将数字频率预畸变转换为模拟频率 Omega = tan(w/2)
# 2. 计算高通滤波器的选择性因子 k = Omega_stop_analog / Omega_pass_analog (注意高通反转)
#    对于高通设计，我们将其映射到低通原型，变换比为 Omega_p / Omega
#    因此选择性因子 k_lp = (tan(wp*pi/2) / tan(ws*pi/2)) 的倒数? 
#    标准低通原型的阻带边缘 Omega_s' = tan(wp*pi/2) / tan(ws*pi/2)

# 数字角频率
w_p_rad = wp * pi
w_s_rad = ws * pi

# 预畸变 (Pre-warping)
Omega_p = tan(w_p_rad / 2)
Omega_s = tan(w_s_rad / 2)

# 计算低通原型的归一化阻带频率 lambda_s
# 对于高通滤波器，通带映射到 1，阻带映射到 Omega_p / Omega_s
lambda_s = Omega_p / Omega_s

println("\n中间计算结果:")
@printf("  模拟通带频率 Ωp: %.4f\n", Omega_p)
@printf("  模拟阻带频率 Ωs: %.4f\n", Omega_s)
@printf("  低通原型阻带比 λs: %.4f\n", lambda_s)

# 切比雪夫阶数公式
# N >= acosh( sqrt((10^(0.1*Rs) - 1) / (10^(0.1*Rp) - 1)) ) / acosh(lambda_s)
numerator = acosh(sqrt((10^(0.1 * Rs) - 1) / (10^(0.1 * Rp) - 1)))
denominator = acosh(lambda_s)
N_exact = numerator / denominator
N = ceil(Int, N_exact)

println("  计算出的最小阶数 N: $N (精确值: $(round(N_exact, digits=3)))")

# ==========================================
# 3. 设计滤波器
# ==========================================
# 使用 DSP.jl 进行设计
# 注意：DSP.jl 的 Highpass(Wn) 中的 Wn 是通带截止频率
# Chebyshev1(N, ripple) 中的 ripple 是通带波纹

response_type = Highpass(wp)
design_method = Chebyshev1(N, Rp)

try
    global filter_obj = digitalfilter(response_type, design_method)
    println("\n滤波器设计成功！")
catch e
    println("\n滤波器设计失败: ", e)
end

# ==========================================
# 4. 计算频率响应并绘图
# ==========================================
# 定义频率向量 (0 到 π)
w_range = range(0, stop=pi, length=512)
h_resp = freqz(filter_obj, w_range)

# 计算幅度 (dB)
mag_db = 20 * log10.(abs.(h_resp))
frequencies_normalized = w_range / pi

println("正在绘制增益响应...")

# 绘制幅频响应
plot(frequencies_normalized, mag_db, "b-", linewidth=2, label="Gain Response")
hold("on")

# 绘制指标限制框
# 1. 通带区域 (0.75 到 1.0): 增益应在 -0.5dB 到 0dB 之间
#    这里只画下限 -0.5dB
plot([wp, 1.0], [-Rp, -Rp], "g--", linewidth=2, label="Passband Ripple Limit (-0.5dB)")

# 2. 阻带区域 (0 到 0.35): 增益应小于 -43dB
plot([0, ws], [-Rs, -Rs], "r--", linewidth=2, label="Stopband Attenuation Limit (-43dB)")

# 绘制截止频率标记
plot([wp, wp], [-60, 5], "k:", label="Passband Edge (0.75\\pi)")
plot([ws, ws], [-60, 5], "k:", label="Stopband Edge (0.35\\pi)")

# 设置图形属性
title("Chebyshev Type I Highpass Filter Response (N=$N)")
xlabel("Normalized Frequency (x \\pi rad/sample)")
ylabel("Magnitude (dB)")
xlim(0, 1)
ylim(-80, 5) # 设置合适的 Y 轴范围以显示衰减
grid("on")
legend()

println("绘图完成。")