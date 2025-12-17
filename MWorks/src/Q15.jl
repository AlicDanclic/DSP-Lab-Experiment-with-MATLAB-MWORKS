using TyPlot
using DSP # Julia的信号处理包，包含filt函数

# 1. 定义系统参数
# H(z) = (1 - 0.2z^-1 + 0.5z^-2) / (1 + 3.2z^-1 + 1.5z^-2 - 0.8z^-3 + 1.4z^-4)
# 分子系数 b (对应 z^-0, z^-1, z^-2 ...)
b = [1.0, -0.2, 0.5]

# 分母系数 a (对应 z^-0, z^-1, z^-2, z^-3, z^-4)
a = [1.0, 3.2, 1.5, -0.8, 1.4]

# 2. 构造输入信号：单位脉冲序列 delta[n]
N = 20           # 样本数量
x = zeros(N)     # 初始化为全0
x[1] = 1.0       # 在 n=0 处 (索引1) 设置为 1

# 3. 计算冲激响应 h[n]
# 使用 filt 函数求解差分方程，等同于 MATLAB 的 filter(b, a, x)
h = filt(b, a, x)

# 打印数值结果
println("前 20 个样本值 h[n]:")
for i in 1:N
    # 打印格式：n=索引值, 样本值保留4位小数
    println("n=$(i-1): $(round(h[i], digits=4))")
end

# 4. 使用 TyPlot 绘图
# 使用 stem 函数绘制离散序列图 (如果 TyPlot 支持 stem)
# 如果环境不支持 stem，可以使用 plot(0:N-1, h, "o-")
figure("Impulse Response") # 创建图形窗口
stem(0:N-1, h, "b-o")      # 绘制火柴杆图 (蓝色，圆点)

title("系统冲激响应 h[n] (前20个点)")
xlabel("采样点 n")
ylabel("幅度")
grid(true)                 # 显示网格