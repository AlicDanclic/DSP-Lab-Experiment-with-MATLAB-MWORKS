using TyPlot
using FFTW
using DSP

# 1. 定义系统参数
# H(e^jw) = B(e^jw) / A(e^jw)
# 分子系数 (对应 e^0, e^-jw, e^-2jw...)
b = [1.0, -0.2, 0.5, 2.0, -0.6]
# 分母系数
a = [1.0, 3.2, 1.5, -0.8, 1.4]

# 2. 计算频率响应
# 设置 FFT 点数（分辨率），点数越多曲线越平滑
N = 1024

# 对系数向量进行补零，使其长度达到 N
# 这是为了利用 FFT 计算离散时间傅里叶变换 (DTFT) 的近似样本
b_padded = [b; zeros(N - length(b))]
a_padded = [a; zeros(N - length(a))]

# 分别计算分子和分母的 FFT
# FFT 计算结果对应频率范围 0 到 2pi (单位圆一圈)
H_num = fft(b_padded)
H_den = fft(a_padded)

# 系统频率响应 H(w) = Num(w) / Den(w)
H_w = H_num ./ H_den

# 3. 提取 0 到 pi 的部分 (通常只需分析前半部分)
half_idx = 1:(div(N, 2) + 1)
H_half = H_w[half_idx]

# 定义对应的频率轴 (0 到 pi)
w_axis = range(0, pi, length=length(half_idx))

# 4. 计算幅度响应和相位响应
# 幅度响应 (dB)
mag_response = 20 * log10.(abs.(H_half))

# 相位响应 (弧度)
phase_response = angle.(H_half)

# 手动实现相位解卷绕 (unwrap) 以获得连续曲线
# 替代 DSP.unwrap 以避免环境依赖和作用域警告
function manual_unwrap(p)
    n = length(p)
    result = copy(p)
    correction = 0.0
    for i in 2:n
        diff = p[i] - p[i-1]
        # 如果相邻点相位跳变超过 pi，则调整 2*pi
        if diff > pi
            correction -= 2 * pi
        elseif diff < -pi
            correction += 2 * pi
        end
        result[i] += correction
    end
    return result
end

phase_response = manual_unwrap(phase_response)

# 5. 绘图
figure("System Frequency Response Analysis")

# 子图 1: 幅频特性
subplot(2, 1, 1)
plot(w_axis ./ pi, mag_response)
title("幅频特性 (Magnitude Response)")
xlabel("归一化频率 (\\times \\pi rad/sample)")
ylabel("幅度 (dB)")
grid("on")

# 子图 2: 相频特性
subplot(2, 1, 2)
plot(w_axis ./ pi, phase_response)
title("相频特性 (Phase Response)")
xlabel("归一化频率 (\\times \\pi rad/sample)")
ylabel("相位 (radians)")
grid("on")

# 调整布局
# tight_layout()