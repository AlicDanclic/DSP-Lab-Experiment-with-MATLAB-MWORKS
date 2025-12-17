using TyPlot
using DSP # 包含 remez 函数和 freqz 函数

# 1. 定义设计指标
Fs = 5000.0       # 采样频率
fp = 1500.0       # 通带截止频率
fs = 1800.0       # 阻带截止频率
dp = 0.015        # 通带波纹
ds = 0.021        # 阻带波纹

# 2. 估算滤波器阶数 (利用 Herrmann 公式)
# 归一化过渡带宽度
Delta_F = (fs - fp) / Fs
log_dp = log10(dp)
log_ds = log10(ds)

# Herrmann 公式系数
a_dp = 0.005309 * (log_dp^2) + 0.07114 * log_dp - 0.4761
g_dp = -0.00266 * (log_dp^2) - 0.5941 * log_dp - 0.4278
D_inf = a_dp * log_ds + g_dp
f_corr = 11.01217 + 0.51244 * (log_dp - log_ds)
N_exact = (D_inf - f_corr * (Delta_F^2)) / Delta_F

# 取整得到滤波器长度 N
N = ceil(Int, N_exact)
# 确保 N 为奇数（对于等波纹低通滤波器，通常建议奇数长度/偶数阶，以便获得更好的特性）
if N % 2 == 0
    N += 1
end

println("估算的滤波器长度 N: $N")

# 3. 使用 Remez 算法设计滤波器
# 权重计算: 权重与波纹成反比
weights = [1.0/dp, 1.0/ds]

# 定义频带: 必须是扁平向量 [start1, stop1, start2, stop2]
# 之前的写法 bands = [(0.0, fp), (fs, Fs/2)] 是错误的
bands = [0.0, fp, fs, Fs/2]

# 定义期望幅度: [1.0, 0.0] (通带为1，阻带为0)
desired = [1.0, 0.0]

# 调用 remez 函数
# 注意: remez(n_taps, bands, desired; weight=..., Hz=...)
h = remez(N, bands, desired, weight=weights, Hz=Fs)

# 4. 计算频率响应
# 使用 freqz 计算复数频率响应
# freqz(b, a, n_points) -> (H, w) or freqz(Filter, range, Fs)
# 这里手动构造频率向量进行计算，或者使用 DSP 的 freqz
f_range = range(0, Fs/2, length=1024)
H_resp = freqz(PolynomialRatio(h, [1.0]), f_range, Fs)

# 计算幅度 (dB)
mag_db = 20 .* log10.(abs.(H_resp))

# 5. 绘图
figure("Remez Filter Design")

# 子图1: 幅频响应 (dB)
subplot(2, 1, 1)
plot(f_range, mag_db, "b-", linewidth=1.5)
title("FIR 低通滤波器幅频响应 (Remez 算法, N=$N)")
ylabel("幅度 (dB)")
xlabel("频率 (Hz)")
grid(true)

# 添加规格限制线以便观察
hold("on")
# 通带下限 (20log10(1-dp)) 和上限 (20log10(1+dp))
pass_lower = 20 * log10(1 - dp)
pass_upper = 20 * log10(1 + dp)
# 阻带上限 (20log10(ds))
stop_limit = 20 * log10(ds)

plot([0, fp], [pass_lower, pass_lower], "g--")
plot([0, fp], [pass_upper, pass_upper], "g--")
plot([fs, Fs/2], [stop_limit, stop_limit], "r--")
text(fp, pass_lower-5, "通带波纹限制")
text(fs, stop_limit+5, "阻带衰减限制")

# 限制 Y 轴范围以便更好地观察细节
ylim(-80, 5)

# 子图2: 冲激响应 h[n]
subplot(2, 1, 2)
stem(0:N-1, h, "k-o", markersize=4)
title("滤波器冲激响应 h[n]")
xlabel("n")
ylabel("幅度")
grid(true)