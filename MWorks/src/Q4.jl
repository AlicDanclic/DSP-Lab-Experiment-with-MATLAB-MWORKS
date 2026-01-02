using DSP
using TyPlot
using Printf

println("=== 切比雪夫 I 型 IIR 高通滤波器设计 ===")

# ==========================================
# 1. 定义性能指标
# ==========================================
Rp_dB = 0.5         # 通带波纹 (dB)
Rs_dB = 43.0        # 最小阻带衰减 (dB)
wp = 0.75           # 通带边缘频率 (归一化, 1.0 = π)
ws = 0.35           # 阻带边缘频率 (归一化, 1.0 = π)

# 将 dB 指标转换为线性幅度指标以便后续绘图和验证
# 通带最小幅度: 10^(-0.5/20)
mag_pass_min = 10^(-Rp_dB / 20)
# 阻带最大幅度: 10^(-43/20)
mag_stop_max = 10^(-Rs_dB / 20)

println("指标参数:")
println("  通带波纹: $Rp_dB dB (线性幅度 > $(round(mag_pass_min, digits=4)))")
println("  阻带衰减: $Rs_dB dB (线性幅度 < $(round(mag_stop_max, digits=5)))")
println("  通带频率: $(wp)π")
println("  阻带频率: $(ws)π")

# ==========================================
# 2. 计算最小阶数 N
# ==========================================
# 使用双线性变换的模拟原型法计算阶数
# 1. 预畸变: 将数字频率映射到模拟频率 Omega = tan(w/2)
Omega_p = tan(wp * pi / 2)
Omega_s = tan(ws * pi / 2)

# 2. 计算高通滤波器的选择性因子 (映射到低通原型)
# 对于高通 HP -> 低通 LP 变换: λ = Ωp / Ω
# 阻带边缘对应的低通原型频率 λs = Ωp / Ωs
lambda_s = Omega_p / Omega_s

# 3. 切比雪夫阶数公式
numerator = acosh(sqrt((10^(0.1 * Rs_dB) - 1) / (10^(0.1 * Rp_dB) - 1)))
denominator = acosh(lambda_s)
N = ceil(Int, numerator / denominator)

println("\n计算结果:")
println("  模拟频率比 λs: $(round(lambda_s, digits=4))")
println("  所需最小阶数 N: $N")

# ==========================================
# 3. 设计滤波器
# ==========================================
# Highpass(Wn) 中的 Wn 是归一化截止频率
# Chebyshev1(N, ripple) 中的 ripple 是 dB
response_type = Highpass(wp)
design_method = Chebyshev1(N, Rp_dB)

filter_obj = digitalfilter(response_type, design_method)
println("  滤波器设计完成。")

# ==========================================
# 4. 计算频率响应并绘图 (线性幅度)
# ==========================================
# 定义频率向量
w_range = range(0, stop=pi, length=1024)
h_resp = freqz(filter_obj, w_range)

# 【优化】计算线性幅度 (Linear Magnitude)
mag = abs.(h_resp)
# 【优化】转换为数组确保 TyPlot 兼容性
frequencies_normalized = collect(w_range ./ pi)

println("正在绘制线性增益响应...")

figure("Chebyshev Type I Highpass Filter", figsize=(10, 8))

# 1. 绘制幅频响应曲线
plot(frequencies_normalized, mag, "b", linewidth=2, label="幅频响应 |H(e^jw)|")
hold("on")

# 2. 绘制指标限制线 (线性坐标)

# (A) 通带区域 (0.75 - 1.0)
# 上限: 1.0 (0dB)
plot([wp, 1.0], [1.0, 1.0], "g--", label="通带上限 (1.0)")
# 下限: 10^(-0.5/20) ≈ 0.944
plot([wp, 1.0], [mag_pass_min, mag_pass_min], "g--", label="通带下限 (-0.5dB)")

# (B) 阻带区域 (0 - 0.35)
# 上限: 10^(-43/20) ≈ 0.007
plot([0, ws], [mag_stop_max, mag_stop_max], "r--", label="阻带上限 (-43dB)")

# (C) 截止频率垂直标记
plot([wp, wp], [0, 1.1], "k:", label="通带截止 (0.75π)")
plot([ws, ws], [0, 1.1], "k:", label="阻带截止 (0.35π)")

# 设置图形属性
title("切比雪夫 I 型高通滤波器 (N=$N) - 线性幅度响应")
xlabel("归一化频率 (×π rad/sample)")
ylabel("幅度 (归一化 0-1)")
xlim(0, 1)
ylim(0, 1.1) # 线性坐标下，稍微多出一点空间看清楚 1.0

grid(true)
legend()

println("绘图完成。")