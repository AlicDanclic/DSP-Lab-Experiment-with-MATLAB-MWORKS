using TyPlot
using FFTW

# 1. 定义滤波器参数
# 6 点滑动平均滤波器
# 系统函数 H(z) = (1/6) * (1 + z^-1 + ... + z^-5)
# 滤波器系数 (单位脉冲响应 h[n])
M = 6
h = ones(M) ./ M

# 2. 计算频率响应
# 使用 FFT 计算离散时间傅里叶变换 (DTFT) 的近似
# 补零到较长的长度以获得平滑的曲线
nfft = 1024
# 补零
h_padded = [h; zeros(nfft - length(h))]

# 计算 FFT
H = fft(h_padded)

# 3. 提取 0 到 pi 的部分 (单边频谱)
# 因为 h[n] 是实数，频谱是共轭对称的，只需分析前半部分
half_len = div(nfft, 2) + 1
H_half = H[1:half_len]

# 定义归一化频率轴 (0 到 1, 对应 0 到 pi)
w_axis = range(0, 1, length=half_len)

# 4. 计算幅频特性和相频特性
# 幅度 (dB)
mag_response = 20 * log10.(abs.(H_half))

# 相位 (弧度)
phase_response = angle.(H_half)

# 手动实现相位解卷绕 (Unwrap)
# 滑动平均滤波器是线性相位滤波器，理论上相位是一条直线
# 但 angle() 函数的结果限制在 [-pi, pi]，会导致跳变
function manual_unwrap(p)
    n = length(p)
    result = copy(p)
    correction = 0.0
    for i in 2:n
        diff = p[i] - p[i-1]
        if diff > pi
            correction -= 2 * pi
        elseif diff < -pi
            correction += 2 * pi
        end
        result[i] += correction
    end
    return result
end

phase_unwrapped = manual_unwrap(phase_response)

# 5. 绘图
figure("6-Point Moving Average Frequency Response")

# 子图 1: 幅频特性
subplot(2, 1, 1)
plot(w_axis, mag_response)
title("6点滑动平均滤波器 - 幅频特性")
xlabel("归一化频率 (\\times \\pi rad/sample)")
ylabel("幅度 (dB)")
grid("on")
# 限制 y 轴范围以避免 -Inf 对显示的影响 (在零点处)
ylim([-50, 5])

# 子图 2: 相频特性
subplot(2, 1, 2)
plot(w_axis, phase_unwrapped)
title("6点滑动平均滤波器 - 相频特性")
xlabel("归一化频率 (\\times \\pi rad/sample)")
ylabel("相位 (radians)")
grid("on")

# 验证线性相位特性：相位应该是一条斜率为 -(M-1)/2 = -2.5 的直线
# 对于 w (0~pi), slope = -2.5. 
# 归一化频率 x = w/pi. 
# 所以 plot(x, phase) 的斜率应该是 -2.5 * pi