using TyPlot
using DSP
using FFTW
# 移除 SpecialFunctions 引用以避免命名冲突

# 1. 滤波器性能指标
# 高通滤波器
ws = 0.45 * pi  # 阻带截止频率
wp = 0.55 * pi  # 通带截止频率
delta = 0.04    # 波纹 (通带和阻带相同)

# 2. 计算参数
# 截止频率 wc (取中心)
wc = (ws + wp) / 2

# 过渡带宽度 dw
dw = abs(wp - ws)

# 阻带衰减 A (dB)
A = -20 * log10(delta)
println("设计参数:")
println("  目标衰减 A = $(round(A, digits=2)) dB")
println("  截止频率 wc = $(wc/pi) * pi")
println("  过渡带 dw = $(dw/pi) * pi")

# 估算阶数 N (使用 Kaiser 经验公式)
# N = (A - 8) / (2.285 * dw)
N_est = (A - 8) / (2.285 * dw)
N = ceil(Int, N_est)
# 保证 N 为偶数 (滤波器长度 M 为奇数)，以保证 Type I FIR (对称，中心不为0)
# 高通滤波器通常需要 Type I (奇数长度 symmetric) 或 Type IV (偶数长度 antisymmetric)
# 这里选择奇数长度 (Type I)
if N % 2 != 0
    N += 1
end
M = N + 1 # 滤波器长度 (Tap 数)
println("  滤波器阶数 N = $N (长度 M = $M)")

# 计算 beta 参数
# 将变量名 beta 修改为 beta_val 避免与函数名冲突
if A > 50
    beta_val = 0.1102 * (A - 8.7)
elseif A >= 21
    beta_val = 0.5842 * (A - 21)^0.4 + 0.07886 * (A - 21)
else
    beta_val = 0.0
end
println("  Kaiser 参数 beta = $(round(beta_val, digits=4))")

# 3. 构造理想高通滤波器 h_d[n]
# 理想高通 = 延迟脉冲 - 理想低通
# h_hp[n] = delta[n - tau] - (wc/pi) * sinc((wc/pi) * (n - tau))
tau = (M - 1) / 2
n = 0:(M-1)

# 注意：Julia 的 sinc(x) 定义为 sin(pi*x)/(pi*x)
# 公式中的 sinc 参数需要归一化
# 理想低通部分系数
h_lp = (wc / pi) .* sinc.((wc / pi) .* (n .- tau))

# 理想高通系数
h_d = -h_lp
# 在中心点 (n = tau) 加上 1 (即 delta 函数)
# 注意索引，Julia 是 1-based，但在 n 向量中我们已经定义了 0 到 M-1
# tau 对应的索引是 tau + 1
center_idx = Int(tau) + 1
h_d[center_idx] += 1.0

# 4. 生成 Kaiser 窗并加窗
w = kaiser(M, beta_val)
h = h_d .* w

# 5. 计算增益响应
nfft = 1024
# 补零
h_padded = [h; zeros(nfft - length(h))]
H = fft(h_padded)
# 取前半部分 (0 到 pi) 进行绘制，这样更直观
half_len = div(nfft, 2) + 1
H_half = H[1:half_len]
# 计算幅度 (dB)
mag_H_db = 20 * log10.(abs.(H_half))

# 频率轴 (0 到 pi)
w_axis = range(0, pi, length=half_len)

# 6. 绘图
figure("High-pass Filter Design using Kaiser Window")

# 子图1: 滤波器系数 h[n]
subplot(2, 1, 1)
stem(n, h)
title("FIR 高通滤波器单位脉冲响应 h[n] (N=$N)")
xlabel("n")
ylabel("幅度")
grid("on")

# 子图2: 增益响应 (dB)
subplot(2, 1, 2)
plot(w_axis ./ pi, mag_H_db)
title("增益响应 (Gain Response) - [0, \\pi]")
xlabel("归一化频率 (\\times \\pi rad/sample)")
ylabel("幅度 (dB)")
xlim([0, 1])     # 显示 0 到 1 (pi)
ylim([-80, 10]) # 限制 Y 轴范围
grid("on")

# 添加指标辅助线
# 阻带截止 (0.45) 和 通带截止 (0.55)
# 阻带上限 (-A dB)
plot([0, 0.45], [-A, -A], "--r", label="阻带指标") 
# 通带下限 (这里近似画在 0dB 附近示意)
plot([0.55, 1], [0, 0], "--g", label="通带指标")
legend()