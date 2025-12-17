import Pkg
# 检查并安装 Polynomials 包，用于多项式运算
if Base.find_package("Polynomials") === nothing
    println("正在安装 Polynomials 包...")
    Pkg.add("Polynomials")
end

using Polynomials
using LinearAlgebra
using Printf

println("=== 系统函数部分分式展开 (ResidueZ) ===")

# ==========================================
# 1. 定义系统参数
# ==========================================
# H(z) = B(z) / A(z)
# 分子系数 (对应 1, z^-1, z^-2 ...)
b = [1.0, -0.2, 0.5]
# 分母系数 (对应 1, z^-1, z^-2, z^-3, z^-4)
a = [1.0, 3.2, 1.5, -0.8, 1.4]

println("分子系数 b: ", b)
println("分母系数 a: ", a)

# ==========================================
# 2. 计算极点 (Poles)
# ==========================================
# 将 H(z) 上下同乘 z^4 (分母最高阶)，转换为 z 的正幂次多项式
# A(z) = 1 + 3.2z^-1 + ... + 1.4z^-4
# -> A_poly(z) = z^4 + 3.2z^3 + 1.5z^2 - 0.8z + 1.4
# 注意：Polynomials 包默认系数顺序为升幂 [a0, a1, a2...]
# 我们的 a 向量是降幂排列 (对应 z^0, z^-1...)，且转换后对应 z^N, z^(N-1)...
# 因此 A_poly 的系数应该是 reverse(a)

den_coeffs = reverse(a) # [1.4, -0.8, 1.5, 3.2, 1.0]
den_poly = Polynomial(den_coeffs)

# 极点 p 是分母的根
# 显式使用 Polynomials.roots 避免与 TyMath.roots 冲突
p = Polynomials.roots(den_poly)

println("\n计算得到的极点 p:")
for (i, val) in enumerate(p)
    @printf("p%d = %.4f %+.4fim\n", i, real(val), imag(val))
end

# ==========================================
# 3. 计算留数 (Residues)
# ==========================================
# 对于 H(z) = B(z)/A(z)，我们需要找到 r_i 使得 H(z) = sum( r_i / (1 - p_i*z^-1) )
# 这等价于对 H(z)/z 进行部分分式展开：H(z)/z = sum( r_i / (z - p_i) )
# r_i = [H(z)/z * (z - p_i)] | z=p_i
#     = Num(p_i) / Den'(p_i)  (对于 H(z)/z 的分子分母)

# 构造 H(z)/z 的分子和分母多项式
# 原 H(z) = (z^4 - 0.2z^3 + 0.5z^2) / (z^4 + ...)  [分子补齐 z^2 以匹配分母阶数]
# H(z)/z = (z^3 - 0.2z^2 + 0.5z) / (z^4 + 3.2z^3 + ...)

# 构造分子系数 (升幂): 0.5z, -0.2z^2, 1.0z^3. 常数项为0
# b 向量是 [1.0, -0.2, 0.5]，对应 1 - 0.2z^-1 + 0.5z^-2
# 乘 z^4 变成 z^4 - 0.2z^3 + 0.5z^2
# 除 z 变成 z^3 - 0.2z^2 + 0.5z
# 系数(升幂): [0.0, 0.5, -0.2, 1.0]
num_coeffs_for_residue = [0.0, 0.5, -0.2, 1.0] 
num_poly_res = Polynomial(num_coeffs_for_residue)

# 分母多项式即为之前的 den_poly (A_poly)
# r_i = Num_res(p_i) / Den_poly'(p_i)

# 显式使用 Polynomials.derivative
den_derivative = Polynomials.derivative(den_poly)

r = ComplexF64[]
for pole in p
    # 计算留数
    res_val = num_poly_res(pole) / den_derivative(pole)
    push!(r, res_val)
end

println("\n计算得到的留数 r:")
for (i, val) in enumerate(r)
    @printf("r%d = %.4f %+.4fim\n", i, real(val), imag(val))
end

# ==========================================
# 4. 生成展开表达式
# ==========================================
println("\n=== 展开后的表达式 ===")
println("H(z) = ")
for i in 1:length(r)
    r_re = real(r[i])
    r_im = imag(r[i])
    p_re = real(p[i])
    p_im = imag(p[i])
    
    # 格式化 r
    r_str = ""
    if abs(r_im) < 1e-6
        r_str = @sprintf("%.4f", r_re)
    else
        r_str = @sprintf("(%.4f %+.4fj)", r_re, r_im)
    end
    
    # 格式化 p
    p_str = ""
    if abs(p_im) < 1e-6
        p_str = @sprintf("1 - %.4fz^-1", p_re)
    else
        p_str = @sprintf("1 - (%.4f %+.4fj)z^-1", p_re, p_im)
    end
    
    # 打印项
    if i == 1
        print("      ", r_str, " / [", p_str, "]")
    else
        print("\n    + ", r_str, " / [", p_str, "]")
    end
end
println("\n")