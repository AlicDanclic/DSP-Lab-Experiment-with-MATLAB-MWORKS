import Pkg
# 检查并安装 TyPlot (如果尚未安装)
if Base.find_package("TyPlot") === nothing
    println("正在安装 TyPlot ...")
    Pkg.add("TyPlot")
end

using TyPlot

println("=== 复指数序列绘制 ===")

# ==========================================
# 1. 定义序列参数
# ==========================================
# 序列范围 n = 0 到 29 (30点)
n = 0:29

# 定义复指数序列 x[n] = 0.2 * e^((0.4 + j0.5)n)
# 在 Julia 中，虚数单位用 im 表示
# 注意使用 .进行广播运算
alpha = 0.4 + 0.5im
x = 0.2 * exp.(alpha .* n)

# ==========================================
# 2. 提取实部和虚部
# ==========================================
x_real = real(x)
x_imag = imag(x)

println("计算完成。准备绘图...")

# ==========================================
# 3. 绘制图形
# ==========================================
# 使用 subplot 将实部和虚部画在同一张图的两个子图中

# --- 子图 1: 实部 ---
subplot(2, 1, 1)
# stem 用于绘制离散序列 (火柴杆图)
# 如果 TyPlot 的 stem 参数与 MATLAB 类似，可以直接传入 x, y
stem(n, x_real, "b") 
title("Real Part of x[n]")
ylabel("Amplitude")
grid("on")

# --- 子图 2: 虚部 ---
subplot(2, 1, 2)
stem(n, x_imag, "r")
title("Imaginary Part of x[n]")
xlabel("n (Time Index)")
ylabel("Amplitude")
grid("on")

println("绘图完成。")