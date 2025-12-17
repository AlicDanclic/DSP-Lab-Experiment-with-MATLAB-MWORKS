# 导入必要的包
using TyPlot
using DSP  # 用于 conv 函数。如果未安装，需先运行 import Pkg; Pkg.add("DSP")

# 1. 定义序列 x1[n] 和 x2[n]
x1 = [2.2, 3.0, -1.5, 4.2, -1.8]
x2 = [0.8, -1.0, 1.6, 0.8]

# 2. 计算卷积 x[n] = x1[n] * x2[n]
# DSP 包中的 conv 函数可以直接计算线性卷积
# 注意：由于 TyMath 和 DSP 都导出了 conv，为了避免冲突，必须显式指定 DSP.conv
xn = DSP.conv(x1, x2)

# 3. 确定时间轴 n
# 假设 x1 和 x2 均从 n=0 开始
# 卷积结果的长度为 L1 + L2 - 1
L = length(xn)
n = 0:(L - 1)

# 4. 打印计算结果
println("卷积结果 x[n]:")
println(xn)

# 5. 使用 TyPlot 绘图
# 初始化图形窗口
figure("Convolution Analysis")

# 使用 stem 绘制离散序列 (火柴杆图)
# 移除了导致报错的 "filled" 参数
stem(n, xn)

# 添加图形装饰
title("序列卷积结果 x[n] = x1[n] * x2[n]")
xlabel("n (样本点)")
ylabel("x[n] 幅度")

# 打开网格
grid("on")

# 保持图形显示 (在某些非交互式环境中可能需要)
# show()