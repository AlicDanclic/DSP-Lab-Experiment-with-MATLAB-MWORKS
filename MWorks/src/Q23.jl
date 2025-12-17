using Polynomials # 用于多项式求根
using Printf

# 1. 定义 FIR 滤波器系数
# H(z) = 2.4 + 3.2z^-1 + 1.5z^-2 + 0.8z^-3 + 1.4z^-4 + 3.6z^-5 + 5.2z^-6
# 对应系数向量 b
b = [2.4, 3.2, 1.5, 0.8, 1.4, 3.6, 5.2]

# 2. 求系统函数的零点
# 构建多项式 P(z) = 2.4*z^6 + ... + 5.2
# 需要反转 b 的顺序传入 Polynomial (因为 Polynomials 系数是从低次到高次)
poly_coeffs = reverse(b)
P = Polynomial(poly_coeffs)
# 使用 Polynomials.roots 求根
zeros_z = Polynomials.roots(P)

# 3. 手动分解为二阶节 (SOS)
# 策略：分离实根和复根，复根按共轭配对，实根两两配对

# 设置容差判断虚部
tol = 1e-5
complex_roots = zeros_z[abs.(imag.(zeros_z)) .> tol]
real_roots = zeros_z[abs.(imag.(zeros_z)) .<= tol]

# 对复根按实部排序，确保共轭对相邻
sort!(complex_roots, by=real)
# 对实根排序
sort!(real_roots)

# 存储配对结果
pairs = []

# 配对复根
for i in 1:2:length(complex_roots)
    if i+1 <= length(complex_roots)
        push!(pairs, (complex_roots[i], complex_roots[i+1]))
    end
end

# 配对实根
for i in 1:2:length(real_roots)
    if i+1 <= length(real_roots)
        push!(pairs, (real_roots[i], real_roots[i+1]))
    end
end

# 4. 打印结果
println("系统函数 H(z) 分解为 3 个二次多项式之积：")
# H(z) = k * H1(z) * H2(z) * H3(z)
# 我们将总增益 k = b[1] = 2.4 分配给第一节
k_total = b[1]

println("形式: H(z) = H1(z) * H2(z) * H3(z)")
println("-"^50)

for (i, pair) in enumerate(pairs)
    r1, r2 = pair
    
    # 计算二阶节系数
    # Section = (1 - r1*z^-1) * (1 - r2*z^-1)
    #         = 1 - (r1+r2)z^-1 + (r1*r2)z^-2
    # 对应的系数为:
    # Constant: 1
    # z^-1: -(r1 + r2)
    # z^-2: r1 * r2
    
    coeff_z1 = -real(r1 + r2)
    coeff_z2 = real(r1 * r2)
    
    # 仅对第一节应用增益 k
    gain = (i == 1) ? k_total : 1.0
    
    c0 = 1.0 * gain
    c1 = coeff_z1 * gain
    c2 = coeff_z2 * gain
    
    # 格式化符号
    sign1 = c1 >= 0 ? "+" : "-"
    sign2 = c2 >= 0 ? "+" : "-"
    
    @printf("H%d(z) = %.4f %s %.4fz^-1 %s %.4fz^-2\n", 
            i, c0, sign1, abs(c1), sign2, abs(c2))
end

println("-"^50)
println("验证提示：各节系数相乘（卷积）应近似等于原系数 [2.4, 3.2, ...]")