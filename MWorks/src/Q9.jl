using TyPlot

# 1. 定义序列范围
# n 从 0 到 23，共 24 个点
n = 0:23

# 2. 计算指数序列 x[n] = 2 * (0.9)^n
# 注意：在 Julia 中，对向量进行幂运算需要使用点运算符 .^
# 这里的 2 .* ... 表示标量与向量的每个元素相乘
x = 2 .* (0.9) .^ n

# 3. 绘制图形
figure("Discrete Exponential Sequence")

# 使用 stem 绘制离散序列 (火柴杆图)
# 语法与 MATLAB 类似
stem(n, x)

# 添加标题和坐标轴标签
# TyPlot 支持 LaTeX 格式的数学公式渲染
title("Exponential Sequence: \$x[n] = 2(0.9)^{n}\$")
xlabel("n (Sample Index)")
ylabel("Amplitude x[n]")

# 开启网格
grid("on")

# 设置 X 轴范围稍微宽一点，以便看清首尾的点
xlim(-1, 24)

println("绘图完成。")