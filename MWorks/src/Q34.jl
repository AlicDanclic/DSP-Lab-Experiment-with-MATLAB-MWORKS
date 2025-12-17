using LinearAlgebra
using TyPlot  # 使用 Syslab 原生绘图库

# ==========================================
#  1. 基础函数：求多项式根
# ==========================================
function get_roots_robust(coeffs)
    # 输入系数向量 [a0, a1, ..., an] 对应 a0 + a1*z^-1 + ...
    n = length(coeffs) - 1
    if abs(coeffs[1]) < 1e-9; return ComplexF64[]; end
    
    # 归一化系数
    c = coeffs ./ coeffs[1]
    
    # 构建伴随矩阵 (Companion Matrix)
    # 特征值即为多项式的根
    Cm = zeros(ComplexF64, n, n)
    for i in 1:n-1; Cm[i+1, i] = 1.0; end
    for i in 1:n; Cm[i, n] = -c[i+1]; end
    
    return eigvals(Cm)
end

# ==========================================
#  2. 定义系统参数
# ==========================================

# 分子系数 (Numerator)
# 1 - 0.2z^-1 + 0.5z^-2 + 2z^-3 - 0.6z^-4
b = [1.0, -0.2, 0.5, 2.0, -0.6]

# 分母系数 (Denominator)
# 1 + 3.2z^-1 + 1.5z^-2 - 0.8z^-3 + 1.4z^-4
a = [1.0, 3.2, 1.5, -0.8, 1.4]

# ==========================================
#  3. 计算零点与极点
# ==========================================

zeros_val = get_roots_robust(b) # 分子根 -> 零点
poles_val = get_roots_robust(a) # 分母根 -> 极点

println("计算完成。")
println("零点 (Zeros): ", round.(zeros_val, digits=4))
println("极点 (Poles): ", round.(poles_val, digits=4))

# ==========================================
#  4. 绘制零极点分布图
# ==========================================

figure("Question 34: Pole-Zero Plot")
hold("on")

# (1) 绘制单位圆 (Unit Circle) - 参考线
theta = range(0, 2pi, length=200)
plot(cos.(theta), sin.(theta), "k--", linewidth=1, label="Unit Circle")

# (2) 绘制零点 (Zeros) - 蓝色圆圈 'o'
plot(real(zeros_val), imag(zeros_val), "bo", markersize=8, label="Zeros")

# (3) 绘制极点 (Poles) - 红色叉号 'x'
plot(real(poles_val), imag(poles_val), "rx", markersize=10, label="Poles")

# (4) 图表设置
title("Pole-Zero Plot of the System")
xlabel("Real Axis")
ylabel("Imaginary Axis")
grid("on")
legend("on")

# 【关键】设置坐标轴比例相等，确保圆画出来是圆的
axis("equal") 

# 如果极点范围很大，适当调整显示范围以便看清
max_val = maximum(abs.([real(poles_val); imag(poles_val); 1.5])) + 0.5
xlim([-max_val, max_val])
ylim([-max_val, max_val])

hold("off")