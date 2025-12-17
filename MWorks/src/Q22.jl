# 导入必要的包
using TyPlot
using FFTW # 用于 fft, fftshift。如果未安装，请运行 import Pkg; Pkg.add("FFTW")

# 1. 定义序列 x[n]
# 范围: 0 <= n <= 15
n = 0:15
# x[n] = cos(pi * n / 2)
# 注意使用 . 进行向量化运算
x = cos.(pi .* n ./ 2)

# 2. 计算 DFT (N点离散傅里叶变换)
# 直接对 16 点序列进行 FFT
X_k = fft(x)
# 取幅度
mag_X_k = abs.(X_k)
# 定义 DFT 的频率索引轴 k
k = 0:length(x)-1

# 3. 计算 DTFT (离散时间傅里叶变换)
# DTFT 是连续频谱。为了在计算机中绘制，我们通常通过对序列进行大量补零
# 然后做 FFT 来获得高密度的频谱样本作为近似。
K = 1024 # 定义高分辨率的点数
# 构造补零后的序列: 原始 x 加上 (K - N) 个零
x_padded = [x; zeros(K - length(x))]

# 计算补零后的 FFT
X_dtft = fft(x_padded)

# 使用 fftshift 将零频率分量移到频谱中心
X_dtft_shifted = fftshift(X_dtft)

# 取幅度
mag_X_dtft = abs.(X_dtft_shifted)

# 定义 DTFT 的频率轴 omega，范围从 -pi 到 pi
w = range(-pi, stop=pi, length=K)

# 4. 绘图
figure("DFT and DTFT Analysis")

# 子图 1: 原始序列 x[n]
subplot(3, 1, 1)
stem(n, x)
title("原始序列 x[n] = cos(\\pi n / 2)")
xlabel("n")
ylabel("x[n]")
grid("on")

# 子图 2: DFT 幅度谱 |X[k]|
subplot(3, 1, 2)
stem(k, mag_X_k)
title("DFT 幅度谱 |X[k]| (N=16)")
xlabel("k (频率索引)")
ylabel("|X[k]|")
grid("on")

# 子图 3: DTFT 幅度谱 |X(e^{j\\omega})|
subplot(3, 1, 3)
plot(w, mag_X_dtft)
title("DTFT 幅度谱 (近似)")
xlabel("\\omega (弧度/样本)")
ylabel("|X(e^{j\\omega})|")
# 设置 x 轴范围为 -pi 到 pi
xlim([-pi, pi])
grid("on")

# 调整布局以防止重叠 (如果支持)
# tight_layout()