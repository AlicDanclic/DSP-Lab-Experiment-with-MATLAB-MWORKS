using TyPlot
using FFTW

# ==========================================
# 0. 核心算法: 手动实现 Kaiser 窗
# ==========================================
# 目的: 避免依赖 DSP.jl 和 SpecialFunctions.jl，确保在 MWorks 中直接运行

function my_besseli0(x)
    s = 1.0; term = 1.0; x_half_sq = (x / 2)^2
    for k in 1:25
        term *= x_half_sq / (k^2); s += term
        if term < 1e-12 * s; break; end
    end
    return s
end

function my_kaiser(M::Int, beta::Float64)
    w = zeros(M); alpha = (M - 1) / 2.0; I0_beta = my_besseli0(beta)
    for n in 0:M-1
        val = beta * sqrt(1 - ((n - alpha) / alpha)^2)
        val = max(0.0, abs(val)) 
        w[n+1] = my_besseli0(val) / I0_beta
    end
    return w
end

# ==========================================
# 1. 滤波器设计计算
# ==========================================
# 指标
wp = 0.55 * pi  # 通带截止
ws = 0.45 * pi  # 阻带截止
delta = 0.04    # 波纹

# 中间参数
wc = (wp + ws) / 2
dw = abs(wp - ws)
A = -20 * log10(delta)

# 计算 Beta
if A > 50; beta = 0.1102 * (A - 8.7)
elseif A >= 21; beta = 0.5842 * (A - 21)^0.4 + 0.07886 * (A - 21)
else; beta = 0.0; end

# 计算阶数 N
N_est = (A - 8) / (2.285 * dw)
N = ceil(Int, N_est)
# 修正: 高通滤波器必须是偶数阶 (奇数长度 Type I)
if N % 2 != 0; N += 1; end 
M = N + 1

println("---------------------------")
println("设计结果:")
println("  阻带衰减 A = $(round(A, digits=2)) dB")
println("  滤波器阶数 N = $N (长度 M = $M)")
println("  Kaiser Beta = $(round(beta, digits=4))")
println("---------------------------")

# 计算系数 h[n]
alpha_val = N / 2
n_seq = 0:N
h_ideal = zeros(M)
for i in 1:M
    val = n_seq[i]
    if val == alpha_val
        h_ideal[i] = 1.0 - (wc / pi)
    else
        m_val = val - alpha_val
        h_ideal[i] = -sin(wc * m_val) / (pi * m_val)
    end
end
# 加窗
h = h_ideal .* my_kaiser(M, beta)

# ==========================================
# 2. 频率响应计算
# ==========================================
K = 2048
H = fft([h; zeros(K - M)])
H_shifted = fftshift(H)

# 频率轴 (转换为数组以确保绘图兼容性)
freqs = collect(range(-1, 1, length=K))

# 幅度响应 (dB)
mag_H = abs.(H_shifted)
gain_dB = 20 .* log10.(mag_H .+ 1e-12) 

# ==========================================
# 3. 绘图结果
# ==========================================
figure("High-pass Filter Design", figsize=(10, 8))

# 子图 1: 冲激响应
subplot(2, 1, 1)
stem(0:N, h, "b-o", label="h[n]")
title("FIR 高通滤波器冲激响应 h[n] (N=$N)")
xlabel("n")
ylabel("幅度")
grid(true)
legend()

# 子图 2: 增益响应
subplot(2, 1, 2)

# 1. 先画指标辅助线 (放在底层)
# 阻带/通带边界
plot([ws/pi, ws/pi], [-100, 20], "k--", label="阻带/通带边界")
plot([wp/pi, wp/pi], [-100, 20], "k--")
# 最小衰减线
plot([0, 1], [-A, -A], "k:", label="最小衰减指标")

# 2. 画增益曲线 (红色实线)
plot(freqs, gain_dB, "r", label="增益响应")

title("增益响应 (dB)")
xlabel("归一化频率 (ω / π)")
ylabel("幅度 (dB)")
xlim(0, 1)      # 只显示 0 到 pi
ylim(-100, 10)  # Y 轴范围: 显示从 -100dB 到 10dB

grid(true)
legend()