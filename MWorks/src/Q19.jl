using TyPlot
using Polynomials # 用于求多项式的根

# 1. 定义系统参数
# H(z) = (1 - 0.2z^-1 + 0.5z^-2) / (1 + 3.2z^-1 + 1.5z^-2 + 0.8z^-3 + 1.4z^-4)

# 将分子分母同时乘以 z^4 转换为 z 的正幂次多项式，以便求根
# Numerator(z)   = z^4 - 0.2z^3 + 0.5z^2 = 0*z^0 + 0*z^1 + 0.5*z^2 - 0.2*z^3 + 1*z^4
# Denominator(z) = z^4 + 3.2z^3 + 1.5z^2 + 0.8z + 1.4 = 1.4*z^0 + 0.8*z^1 + 1.5*z^2 + 3.2*z^3 + 1*z^4

# Polynomials 包的系数顺序通常为升幂排列 [a0, a1, ..., an] 对应 a0 + a1*x + ...
# 分子系数 (升幂)
b_coeffs_asc = [0.0, 0.0, 0.5, -0.2, 1.0]
# 分母系数 (升幂)
a_coeffs_asc = [1.4, 0.8, 1.5, 3.2, 1.0]

# 2. 计算零点和极点
# 构造多项式对象
num_poly = Polynomial(b_coeffs_asc)
den_poly = Polynomial(a_coeffs_asc)

# 求根
zeros_z = Polynomials.roots(num_poly)
poles_z = Polynomials.roots(den_poly)

# 3. 稳定性判断
# 计算极点模长
poles_abs = abs.(poles_z)

# 判断是否所有极点都在单位圆内
is_stable = all(poles_abs .< 1.0)

# 输出结果
println("---------- 系统稳定性分析 ----------")
println("系统极点 (Poles):")
for p in poles_z
    println("  $p (模长: $(round(abs(p), digits=4)))")
end

println("\n最大极点模长: $(round(maximum(poles_abs), digits=4))")

if is_stable
    println("结论: 系统是稳定的 (所有极点在单位圆内)。")
else
    println("结论: 系统是不稳定的 (存在极点在单位圆外)。")
end
println("------------------------------------")

# 4. 绘制零极点图 (Pole-Zero Plot)
figure("Pole-Zero Plot")
hold("on") # 修正：使用 "on" 字符串而不是布尔值

# 绘制单位圆
theta = 0:0.01:2*pi
unit_circle_x = cos.(theta)
unit_circle_y = sin.(theta)
plot(unit_circle_x, unit_circle_y, "k--", label="单位圆")

# 绘制零点 (圆圈 o)
plot(real(zeros_z), imag(zeros_z), "bo", markersize=8, label="零点 (Zeros)")

# 绘制极点 (叉号 x)
plot(real(poles_z), imag(poles_z), "rx", markersize=8, label="极点 (Poles)")

title("系统零极点分布图")
xlabel("实部 (Real)")
ylabel("虚部 (Imaginary)")
legend()
grid(true)
axis("equal") # 保持横纵坐标比例一致，使圆看起来是圆的