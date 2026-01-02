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

# 6. 绘制规格示意图（更标准的阶梯/边界显示）
figure("Filter Specifications")

# 一次 plot 画多段线（不使用 hold，避免版本差异）
plot([0, fp],      [1, 1], "r--",      # 通带
     [fp, fs],     [1, 0], "r--",      # 过渡带（示意）
     [fs, Fs/2],   [0, 0], "r--",      # 阻带
     [fp, fp],     [0, 1], "k:",       # fp 竖线
     [fs, fs],     [0, 1], "k:")       # fs 竖线

title("FIR 低通滤波器设计规格")
xlabel("频率 (Hz)")
ylabel("理想幅度")
grid("on")
axis("tight")
xlim([0, Fs/2])
ylim([-0.1, 1.1])

# 标注位置下移，避免挡标题
text(fp, 0.92, "fp=1500 Hz")
text(fs, 0.08, "fs=1800 Hz")