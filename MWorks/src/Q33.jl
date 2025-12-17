using TyPlot # 使用 Syslab 原生绘图库
using LinearAlgebra

# ==========================================
#  1. 基础函数定义
# ==========================================

# 手动实现 FIR 滤波函数 (卷积)
function my_fir_filter(b, x)
    N = length(x)
    M = length(b)
    y = zeros(Float64, N)
    
    for n in 1:N
        # y[n] = sum(b[k] * x[n-k])
        # 注意 Julia 索引从 1 开始，所以对应公式中的 n-1, n-2...
        val = 0.0
        for k in 1:M
            if n - k + 1 > 0
                val += b[k] * x[n - k + 1]
            end
        end
        y[n] = val
    end
    return y
end

# ==========================================
#  2. 系统参数与信号生成
# ==========================================

# 滤波器系数 (根据差分方程)
# y[n] = -6.76195x[n] + 13.456335x[n-1] - 6.76195x[n-2]
b = [-6.76195, 13.456335, -6.76195]

# 时间轴 (取足够长以观察包络)
# w=0.1 对应的周期是 2pi/0.1 ≈ 62.8 点
# 我们取 0 到 150 点，大约包含 2.5 个低频周期
n = 0:200 

# 生成输入信号
# x[n] = [cos(0.1n) + cos(0.4n)] * u[n]
x = [cos(0.1 * i) + cos(0.4 * i) for i in n]

# ==========================================
#  3. 执行滤波
# ==========================================

y = my_fir_filter(b, x)

# ==========================================
#  4. 绘图 (输入与输出的包络)
# ==========================================

figure("Question 33: Envelope Analysis")

# --- 子图 1: 输入序列 x[n] ---
subplot(2, 1, 1)
# 为了展示"包络"效果，我们同时画出离散杆图和连续连线
stem(n, x, "b-", markersize=2, label="Discrete Samples") # 离散点
hold("on")
plot(n, x, "b--", linewidth=1, alpha=0.5, label="Envelope") # 连线表示包络趋势
title("Input Signal x[n]: cos(0.1n) + cos(0.4n)")
xlabel("n")
ylabel("Amplitude")
grid("on")
legend("on")

# --- 子图 2: 输出序列 y[n] ---
subplot(2, 1, 2)
stem(n, y, "r-", markersize=2, label="Discrete Samples")
plot(n, y, "r--", linewidth=1, alpha=0.5, label="Envelope")
title("Output Signal y[n] (0.1 rad/s component filtered out)")
xlabel("n")
ylabel("Amplitude")
grid("on")
legend("on")

hold("off")

# ==========================================
#  5. 验证分析 (打印到控制台)
# ==========================================
println("================ 分析结果 ================")
println("输入信号包含频率: 0.1 rad 和 0.4 rad")
println("观察输出图形:")
println("1. 初始阶段有短暂的瞬态响应。")
println("2. 稳态后，输出变成单一频率的正弦波。")
println("   (因为滤波器是一个陷波器，滤除了 0.1 rad 的分量)")