using TyPlot
using DSP
using Random

# 1. 生成信号和噪声
# 设定序列长度
N = 60
n = 0:(N-1)

# 原始序列 s[n] = 3 * n * (0.8)^n
s = 3 .* n .* (0.8).^n

# 生成加性噪声 d[n]
# "幅度 0.6" 这里理解为噪声在 [-0.6, 0.6] 之间均匀分布
Random.seed!(123) # 固定随机种子以保证结果可复现
d = 0.6 .* (2 .* rand(N) .- 1)

# 受干扰的序列 x[n]
x = s + d

# 2. 定义滤波器
# 5 点滑动平均滤波器
# h[n] = [1/5, 1/5, 1/5, 1/5, 1/5]
M = 5
h = ones(M) ./ M

# 3. 使用 fftfilt 进行滤波
# fftfilt 函数利用重叠相加法(Overlap-Add)高效计算长序列的卷积
y = fftfilt(h, x)

# 4. 绘图
figure("Overlap-Add Filtering using fftfilt")

# 子图 1: 原始未受干扰序列 s[n]
subplot(3, 1, 1)
stem(n, s)
title("原始序列 s[n] = 3n(0.8)^n")
xlabel("n")
ylabel("幅度")
grid("on")

# 子图 2: 受干扰的序列 x[n] (为了对比展示)
subplot(3, 1, 2)
stem(n, x)
title("受干扰序列 x[n] = s[n] + d[n]")
xlabel("n")
ylabel("幅度")
grid("on")

# 子图 3: 滤波器输出的有噪序列 y[n]
subplot(3, 1, 3)
stem(0:(length(y)-1), y)
title("5点滑动平均滤波器输出 (利用 fftfilt)")
xlabel("n")
ylabel("幅度")
grid("on")

# 调整布局
# tight_layout()