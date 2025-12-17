using TyPlot       # Syslab 原生绘图
using LinearAlgebra
using Printf

# ==========================================
#  1. 辅助函数：从根构建多项式系数
# ==========================================
# 相当于 MATLAB 的 poly() 函数
function poly_from_roots(roots_val)
    # 初始化多项式为 [1] (对应最高次项系数)
    coeffs = [1.0 + 0.0im]
    
    for r in roots_val
        # 多项式卷积: (coeff...)*z - (coeff...)*r
        # 原系数左移一位 (乘z) + 原系数乘-r
        new_len = length(coeffs) + 1
        new_coeffs = zeros(ComplexF64, new_len)
        
        new_coeffs[1:end-1] += coeffs       # z 的高次项
        new_coeffs[2:end]   -= coeffs .* r   # 常数项/低次项
        
        coeffs = new_coeffs
    end
    
    # 由于共轭根的存在，理论上虚部应抵消为0，取实部即可
    return real.(coeffs)
end

# ==========================================
#  2. 定义零点和极点
# ==========================================
# 零点 Zeros
z_roots = [
    2.2, 
    -1.0 + 1.0im, 
    -1.0 - 1.0im, 
    1.4
]

# 极点 Poles
p_roots = [
    3.7 + 2.0im, 
    3.7 - 2.0im, 
    -2.1 - 1.0im, 
    -2.1 + 1.0im
]

# ==========================================
#  3. 计算系统函数系数
# ==========================================
b = poly_from_roots(z_roots) # 分子系数
a = poly_from_roots(p_roots) # 分母系数

println("========================================")
println("          第 38 题计算结果")
println("========================================")

println("系统函数 H(z) 的系数 (标准形式):")
println("H(z) = (b0 + b1*z^-1 + ... ) / (a0 + a1*z^-1 + ...)")
println("----------------------------------------")

println("[分子系数 b]:")
println(round.(b, digits=4))

println("\n[分母系数 a]:")
println(round.(a, digits=4))

println("\n数学表达式 H(z):")
print("H(z) = ( 1")
@printf(" %.4f z^-1 + %.4f z^-2 + %.4f z^-3 + %.4f z^-4 )\n", b[2], b[3], b[4], b[5])
print("       -------------------------------------------------------------------------\n")
print("       ( 1")
@printf(" %.4f z^-1 + %.4f z^-2 + %.4f z^-3 + %.4f z^-4 )\n", a[2], a[3], a[4], a[5])


# ==========================================
#  4. 绘制零极点图 (验证用)
# ==========================================
figure("Question 38: Pole-Zero Plot")
hold("on")

# 画单位圆
theta = range(0, 2pi, length=200)
plot(cos.(theta), sin.(theta), "k--", linewidth=1.5, label="Unit Circle")

# 画零点
plot(real(z_roots), imag(z_roots), "bo", markersize=9, label="Zeros")
# 画极点
plot(real(p_roots), imag(p_roots), "rx", markersize=10, label="Poles")

title("Pole-Zero Plot")
xlabel("Real Axis")
ylabel("Imaginary Axis")
grid("on")
legend("on")
axis("equal")

# 调整视野以包含所有点
max_val = maximum(abs.([z_roots; p_roots])) * 1.2
xlim([-max_val, max_val])
ylim([-max_val, max_val])

hold("off")