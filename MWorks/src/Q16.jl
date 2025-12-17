using TyPlot

# 1. 定义滤波器指标
Fs = 5000.0       # 抽样频率 (Hz)
fp = 1500.0       # 通带截止频率 (Hz)
fs = 1800.0       # 阻带截止频率 (Hz)
dp = 0.015        # 通带波纹 (delta_p)
ds = 0.021        # 阻带波纹 (delta_s)

# 2. 预处理参数
# 计算归一化过渡带宽度 Delta F
# 注意：Herrmann 公式中的频率通常归一化为 Fs=1
Delta_F = (fs - fp) / Fs

# 计算对数值
log_dp = log10(dp)
log_ds = log10(ds)

# 3. 计算 Hermann 公式系数
# 系数定义 (参考 Herrmann, Schuessler, Dehner, 1973)
# D_inf = a(dp) * log10(ds) + g(dp)

# 计算 a(dp)
# a(dp) = 0.005309 * (log_dp)^2 + 0.07114 * log_dp - 0.4761
a_dp = 0.005309 * (log_dp^2) + 0.07114 * log_dp - 0.4761

# 计算 g(dp)
# g(dp) = -0.00266 * (log_dp)^2 - 0.5941 * log_dp - 0.4278
g_dp = -0.00266 * (log_dp^2) - 0.5941 * log_dp - 0.4278

# 计算 D_infinity
D_inf = a_dp * log_ds + g_dp

# 计算修正项 f(dp, ds)
# f(dp, ds) = 11.01217 + 0.51244 * (log_dp - log_ds)
f_correction = 11.01217 + 0.51244 * (log_dp - log_ds)

# 4. 计算滤波器长度 N
# N = (D_inf - f_correction * (Delta_F)^2) / Delta_F
N_exact = (D_inf - f_correction * (Delta_F^2)) / Delta_F

# 滤波器长度通常取向上取整
N = ceil(Int, N_exact)

# 滤波器阶数 M = N - 1
Order = N - 1

# 5. 输出结果
println("---------- Hermann 公式估算结果 ----------")
println("归一化过渡带宽度 ΔF: $(round(Delta_F, digits=4))")
println("D_infinity: $(round(D_inf, digits=4))")
println("计算得到的精确长度 N_exact: $(round(N_exact, digits=4))")
println("估算滤波器长度 N: $N")
println("估算滤波器阶数 M (N-1): $Order")
println("----------------------------------------")

# 6. (可选) 绘制规格示意图
# 这是一个简单的图示，用于直观展示通带和阻带的位置
figure("Filter Specifications")
# 绘制理想的幅频特性轮廓
# 修改变量名 freqs 为 freq_specs 以避免与 DSP.freqs 冲突
freq_specs = [0, fp, fs, Fs/2]
amps = [1, 1, 0, 0]
plot(freq_specs, amps, "r--")
title("FIR 低通滤波器设计规格")
xlabel("频率 (Hz)")
ylabel("理想幅度")
text(fp, 1.05, "fp=1500")
text(fs, 0.1, "fs=1800")
grid(true)