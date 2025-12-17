using DSP
using TyPlot

# ===========================
# 1. 定义设计指标
# ===========================
Fs = 5000.0         # 抽样频率 (Hz)
fp = 1000.0         # 通带截止频率 (Hz)
fs_stop = 1250.0    # 阻带截止频率 (Hz)
Rp = 0.4            # 通带波纹 (dB)
Rs = 45.0           # 最小阻带衰减 (dB)

# ===========================
# 2. 滤波器设计
# ===========================
nyquist = Fs / 2
wp = fp / nyquist
ws = fs_stop / nyquist

println("正在计算滤波器参数...")
# 计算阶数
N, Wn = ellipord(wp, ws, Rp, Rs, domain=:z)
println("滤波器阶数 N: ", N)

# 设计滤波器
# ！！！修改处：显式使用 DSP.Elliptic 以避免命名冲突 ！！！
filt_design = digitalfilter(Lowpass(Wn), DSP.Elliptic(N, Rp, Rs))

# ===========================
# 3. 计算频率响应
# ===========================
nb_points = 1024
freq_range = range(0, nyquist, length=nb_points)

# 计算复数响应
H = freqz(filt_design, freq_range, Fs)

# 转换为 dB
mag_db = 20 * log10.(abs.(H))

# ===========================
# 4. 使用 TyPlot 绘图
# ===========================
println("正在使用 TyPlot 绘图...")

# 创建图形窗口
figure()

# 1. 绘制主响应曲线 (蓝色实线)
plot(freq_range, mag_db, "b", linewidth=1.5)

# 保持图形，以便叠加辅助线
hold("on")

# 2. 绘制通带波纹限制线 (红色虚线)
plot([0, fp], [-Rp, -Rp], "r--", linewidth=1.5)

# 3. 绘制阻带衰减限制线 (绿色虚线)
plot([fs_stop, nyquist], [-Rs, -Rs], "g--", linewidth=1.5)

# 4. 绘制截止频率竖线 (黑色点线，辅助观察)
plot([fp, fp], [-100, 5], "k:")
plot([fs_stop, fs_stop], [-100, 5], "k:")

# 设置图形属性
title("椭圆 IIR 数字低通滤波器增益响应 (TyPlot)")
xlabel("频率 (Hz)")
ylabel("幅度 (dB)")

# 限制 Y 轴范围方便观察
ylim([-100, 5])
xlim([0, nyquist])

# 打开网格
grid("on")

# 添加图例
legend(["滤波器响应", "通带限制 (-0.4dB)", "阻带限制 (-45dB)", "截止频率标记"])

println("绘图完成。")