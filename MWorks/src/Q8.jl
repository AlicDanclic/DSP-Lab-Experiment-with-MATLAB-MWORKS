using DSP
using TyPlot
using Printf

println("=== 巴特沃兹 IIR 带通滤波器设计 ===")

# ==========================================
# 1. 定义性能指标
# ==========================================
# 归一化频率 (1.0 = π)
wp = (0.4, 0.6)      # 通带范围
ws = (0.3, 0.7)      # 阻带范围
Rp_dB = 0.6          # 通带波纹 (dB)
Rs_dB = 35.0         # 最小阻带衰减 (dB)

# 将 dB 指标转换为线性幅度以便绘图
# 通带最小幅度: 10^(-0.6/20)
mag_pass_min = 10^(-Rp_dB / 20)
# 阻带最大幅度: 10^(-35/20)
mag_stop_max = 10^(-Rs_dB / 20)

println("指标参数:")
println("  通带频率: $(wp[1])π 到 $(wp[2])π")
println("  阻带频率: <$(ws[1])π 和 >$(ws[2])π")
println("  通带波纹: $Rp_dB dB (线性幅度 > $(round(mag_pass_min, digits=4)))")
println("  阻带衰减: $Rs_dB dB (线性幅度 < $(round(mag_stop_max, digits=5)))")

# ==========================================
# 2. 计算阶数和截止频率
# ==========================================
# buttord 返回最小阶数 N 和 3dB 截止频率 Wn
N, Wn = buttord(wp, ws, Rp_dB, Rs_dB)

println("\n设计结果:")
println("  滤波器阶数 N: $N")
println("  3dB 截止频率 Wn: $(round.(Wn, digits=4))")

# ==========================================
# 3. 设计滤波器
# ==========================================
responsetype = Bandpass(Wn[1], Wn[2])
designmethod = Butterworth(N)
filter_system = digitalfilter(responsetype, designmethod)

println("  滤波器设计完成。")

# ==========================================
# 4. 计算频率响应并绘图 (线性幅度)
# ==========================================
# 定义频率向量
w_range = range(0, stop=pi, length=1024)
h_resp = freqz(filter_system, w_range)

# 【优化】计算线性幅度
mag = abs.(h_resp)
# 【优化】转换为数组确保 TyPlot 兼容性
frequencies_normalized = collect(w_range ./ pi)

println("正在绘制线性增益响应...")

figure("Butterworth IIR Bandpass Filter", figsize=(10, 8))

# 1. 绘制幅频响应曲线
plot(frequencies_normalized, mag, "b", linewidth=2, label="幅频响应 |H(e^jw)|")
hold("on")

# 2. 绘制指标限制线 (线性坐标)

# (A) 通带区域 (0.4 - 0.6)
# 上限: 1.0
plot([wp[1], wp[2]], [1.0, 1.0], "g--", label="通带上限 (1.0)")
# 下限: 10^(-0.6/20)
plot([wp[1], wp[2]], [mag_pass_min, mag_pass_min], "g--", label="通带下限 (-0.6dB)")

# (B) 阻带区域 1 (0 - 0.3)
plot([0, ws[1]], [mag_stop_max, mag_stop_max], "r--", label="阻带上限 (-35dB)")

# (C) 阻带区域 2 (0.7 - 1.0)
plot([ws[2], 1.0], [mag_stop_max, mag_stop_max], "r--")

# (D) 垂直标记线
plot([wp[1], wp[1]], [0, 1.1], "k:", label="通带边界")
plot([wp[2], wp[2]], [0, 1.1], "k:")
plot([ws[1], ws[1]], [0, 1.1], "k--", alpha=0.5, label="阻带边界")
plot([ws[2], ws[2]], [0, 1.1], "k--", alpha=0.5)

# 设置图形属性
title("巴特沃兹 IIR 带通滤波器 (N=$N) - 线性幅度响应")
xlabel("归一化频率 (×π rad/sample)")
ylabel("幅度 (归一化 0-1)")
xlim(0, 1)
ylim(0, 1.1) 

grid(true)
legend()

println("绘图完成。")