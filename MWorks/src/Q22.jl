using TyPlot
using FFTW

# ==========================================
# 1. 信号定义
# ==========================================
# 题目: x[n] = cos(pi * n / 2), 0 <= n <= 15
N = 16
n = 0:N-1
x = cos.(pi .* n ./ 2)

# ==========================================
# 2. 计算 DFT (离散傅里叶变换)
# ==========================================
# DFT 计算的是离散频点 k = 0, 1, ..., N-1
X_k = fft(x)
mag_X_k = abs.(X_k)
k = 0:N-1

# ==========================================
# 3. 计算 DTFT (离散时间傅里叶变换) - 近似
# ==========================================
# DTFT 是连续频谱 X(e^jw)。
# 为了在计算机中绘制，我们对信号进行大量补零，
# 然后使用 FFT 计算出高密度的频域样本。

# 设置高分辨率点数 (比如 1024 或 2048)
K = 2048 

# 构造补零后的信号: [原始信号, 0, 0, ...]
x_padded = [x; zeros(K - N)]

# 计算 FFT 并移位，使零频率位于中心
X_dtft = fft(x_padded)
X_dtft_shifted = fftshift(X_dtft)
mag_X_dtft = abs.(X_dtft_shifted)

# 生成归一化频率轴 w / pi
# 范围: -1 到 1 (对应 -pi 到 pi)
# 这样画图时，x轴为 0.5 就代表 0.5pi (即 pi/2)
w_normalized = range(-1, 1, length=K)

# ==========================================
# 4. 绘图
# ==========================================
figure("DFT and DTFT Analysis", figsize=(10, 8))

# 子图 1: 时域序列 x[n]
subplot(3, 1, 1)
# 修正: 移除了不兼容的 basefmt 参数
stem(n, x, "b-o", label="x[n]")
title("时域序列 x[n] = cos(\\pi n / 2)")
xlabel("n")
ylabel("幅度")
grid(true)
xlim([-1, 16])

# 子图 2: DFT 幅度谱 |X[k]| (N点)
subplot(3, 1, 2)
# 修正: 移除了不兼容的 basefmt 参数
stem(k, mag_X_k, "r-o", label="|X[k]|")
title("DFT 幅度谱 |X[k]| (N=16)")
xlabel("频率索引 k")
ylabel("幅度")
grid(true)
xlim([0, 15])
# 解释: 对于 cos(pi*n/2)，频率是 1/4 fs。
# k = N * f = 16 * 0.25 = 4。所以 k=4 和 k=12 (16-4) 处有峰值。

# 子图 3: DTFT 幅度谱 |X(e^{j\\omega})| (归一化频率)
subplot(3, 1, 3)
plot(w_normalized, mag_X_dtft, "g", linewidth=2, label="|X(e^{j\\omega})|")
title("DTFT 幅度谱 (频率轴归一化: \\times \\pi rad/sample)")
xlabel("归一化频率 (\\omega / \\pi)") # 这样轴刻度 0.5 就代表 pi/2
ylabel("幅度")
grid(true)
xlim([-1, 1])

# 叠加显示 DFT 的点在 DTFT 上对应的位置 (为了展示 DFT 是 DTFT 的采样)
# DFT 的 k 对应模拟频率 w_k = 2*pi*k/N
# 归一化频率 w_norm = w_k / pi = 2*k/N
# 注意 fftshift 后的范围问题，为了简单展示，这里只画 DTFT 曲线
# 峰值应该出现在 +/- 0.5 处 (因为 pi/2 / pi = 0.5)

# 自动调整布局
# tight_layout()