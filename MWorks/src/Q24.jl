using TyPlot
using DSP
using FFTW

# 1. 滤波器性能指标
wp = 0.35 * pi  # 通带截止频率
ws = 0.45 * pi  # 阻带截止频率
delta = 0.01    # 波纹 (通带和阻带相同)

# 2. 计算 Kaiser 窗参数
# 2.1 计算最小阻带衰减 A (dB)
A = -20 * log10(delta)
println("目标衰减量 A = $A dB")

# 2.2 计算过渡带宽度 dw (rad)
dw = ws - wp
println("过渡带宽度 dw = $(dw/pi) * pi")

# 2.3 估算滤波器阶数 N (经验公式)
# 公式: N >= (A - 8) / (2.285 * dw)
N_float = (A - 8) / (2.285 * dw)
N = ceil(Int, N_float) 
# 确保 N 为偶数 (通常 Type I FIR 滤波器阶数为偶数，长度为奇数)
if N % 2 != 0
    N += 1
end
L = N + 1 # 窗函数长度
println("估算滤波器阶数 N = $N (窗长 L = $L)")

# 2.4 计算形状参数 beta
# 公式:
# if A > 50: beta = 0.1102(A - 8.7)
# if 21 <= A <= 50: beta = 0.5842(A - 21)^0.4 + 0.07886(A - 21)
# if A < 21: beta = 0
if A > 50
    beta = 0.1102 * (A - 8.7)
elseif A >= 21
    beta = 0.5842 * (A - 21)^0.4 + 0.07886 * (A - 21)
else
    beta = 0.0
end
println("Kaiser 窗参数 beta = $beta")

# 3. 生成 Kaiser 窗函数 w[n]
# DSP 包中的 kaiser 函数生成长度为 L 的窗
w = kaiser(L, beta)

# 4. 计算窗函数的频率响应 (增益响应)
# 为了获得光滑的频谱曲线，进行大量补零 (例如补到 1024 点)
nfft = 1024

# 修复：Julia 的 fft 函数第二个参数是变换维度，不是 FFT 长度。
# 若要进行 N 点 FFT（补零），必须手动先对数组进行补零。
w_padded = [w; zeros(nfft - length(w))]

# 对补零后的序列进行 FFT
W = fft(w_padded)
W_shifted = fftshift(W) # 将零频移到中心
mag_W = abs.(W_shifted)

# 归一化幅度 (最大值为 0 dB)
# 注意：对于窗函数本身的频谱，通常看其主瓣宽度和旁瓣衰减
# 这里我们将直流分量归一化到 0dB 以便观察相对衰减
mag_W_db = 20 * log10.(mag_W ./ maximum(mag_W))

# 频率轴 (-pi 到 pi)
freq_axis = range(-pi, pi, length=nfft)

# 5. 绘图
figure("Kaiser Window Design")

# 子图1: 时域波形 w[n]
subplot(2, 1, 1)
n_seq = 0:(L-1)
stem(n_seq, w)
title("Kaiser 窗函数时域波形 (N=$N, \\beta=$(round(beta, digits=2)))")
xlabel("n")
ylabel("幅度")
grid("on")

# 子图2: 频域增益响应 (dB)
subplot(2, 1, 2)
plot(freq_axis ./ pi, mag_W_db)
title("Kaiser 窗函数幅频响应 (增益)")
xlabel("归一化频率 (\\times \\pi rad/sample)")
ylabel("幅度 (dB)")
xlim([-1, 1])
ylim([-100, 10]) # 限制 Y 轴范围以便观察旁瓣
grid("on")

# 标记主瓣和旁瓣衰减的辅助线 (可选)
# 理论旁瓣峰值应该在 -A dB 左右 (相对于归一化后的主瓣其实是衡量相对值，
# 但对于窗函数谱，第一旁瓣通常由 beta 决定)