using Polynomials
using Printf

# 1. 定义系统函数系数
# H(z) = B(z) / A(z)
# 分子系数 b: 2 + 5z^-1 + 1z^-2 - 3z^-3 + 4z^-4 + 6z^-5
b = [2.0, 5.0, 1.0, -3.0, 4.0, 6.0]
# 分母系数 a: 1 + 3z^-1 - 5z^-2 + 2z^-3 - 4z^-4 + 3z^-5
a = [1.0, 3.0, -5.0, 2.0, -4.0, 3.0]

println("原系统函数 H(z):")
println("分子 b: ", b)
println("分母 a: ", a)
println("-"^60)

# 2. 求零点和极点
# 为了求 z 的根，我们需要构建多项式。
# H(z) 的各项是 z^-k，乘以 z^N 后变成关于 z 的多项式。
# 系数向量需要反转（因为 Polynomials 默认系数顺序是 [z^0, z^1, ..., z^N]）
zeros_z = Polynomials.roots(Polynomial(reverse(b)))
poles_p = Polynomials.roots(Polynomial(reverse(a)))

# 3. 零极点配对 (Pairing) 算法
# 目标：将零点和极点分组，每组包含 2 个零点和 2 个极点（不足补0），形成二阶节
function pair_roots_complex(r_list)
    # 分离实根和复根
    tol = 1e-5
    complex_r = r_list[abs.(imag.(r_list)) .> tol]
    real_r = r_list[abs.(imag.(r_list)) .<= tol]
    
    # 简单的排序策略：按实部排序
    sort!(complex_r, by=real)
    # 修改：Julia 无法直接比较复数大小，即使它们近似为实数。
    # 必须显式指定按实部 (by=real) 进行排序。
    sort!(real_r, by=real)
    
    pairs = []
    
    # 配对复根 (共轭对)
    i = 1
    while i < length(complex_r)
        push!(pairs, (complex_r[i], complex_r[i+1]))
        i += 2
    end
    # 如果剩下一个复根 (理论上实系数多项式不会发生，除非数值误差)，暂且放入实根处理
    if i == length(complex_r)
        push!(real_r, complex_r[end])
    end
    
    # 配对实根
    j = 1
    while j < length(real_r)
        push!(pairs, (real_r[j], real_r[j+1]))
        j += 2
    end
    
    # 如果还剩一个实根 (奇数阶情况)
    if j == length(real_r)
        # 配对一个 0 (在 z 平面原点，意味着没有该项)
        # 注意：这里我们只存储根。之后转换系数时，单个根 r 对应 (z-r)
        push!(pairs, (real_r[end], nothing))
    end
    
    return pairs
end

zero_pairs = pair_roots_complex(zeros_z)
pole_pairs = pair_roots_complex(poles_p)

# 确保零点对和极点对数量一致 (不足的补 1，即 (z-0)(z-0) 对应的系数是 1)
# 这里假设分子分母阶数相同或相近
n_sections = max(length(zero_pairs), length(pole_pairs))

# 4. 构建级联结构 (SOS) 并输出
# 计算总增益 k = b[0] / a[0]
k_gain = b[1] / a[1]

println("级联结构分解结果 (H(z) = H1(z) * H2(z) * ...):")
println("注意: 第一个子系统包含了总增益 k = $(k_gain)")
println("-"^60)

for i in 1:n_sections
    # 获取零点对
    z_pair = i <= length(zero_pairs) ? zero_pairs[i] : (0.0, 0.0)
    # 获取极点对
    p_pair = i <= length(pole_pairs) ? pole_pairs[i] : (0.0, 0.0)
    
    # --- 计算分子系数 (Numerator) ---
    # (1 - z1*z^-1)(1 - z2*z^-1) = 1 - (z1+z2)z^-1 + (z1*z2)z^-2
    if z_pair[2] === nothing
        # 单实根情况: (1 - z1*z^-1)
        z1 = z_pair[1]
        b_sec = [1.0, -real(z1), 0.0]
    else
        z1, z2 = z_pair
        b_sec = [1.0, -real(z1 + z2), real(z1 * z2)]
    end
    
    # --- 计算分母系数 (Denominator) ---
    if p_pair[2] === nothing
        p1 = p_pair[1]
        a_sec = [1.0, -real(p1), 0.0]
    else
        p1, p2 = p_pair
        a_sec = [1.0, -real(p1 + p2), real(p1 * p2)]
    end
    
    # --- 分配增益 ---
    if i == 1
        b_sec = b_sec .* k_gain
    end
    
    # --- 格式化输出 ---
    # 辅助函数：格式化项
    function fmt_term(val, power)
        if abs(val) < 1e-6 return "" end # 忽略极小值
        sign_str = val >= 0 ? "+" : "-"
        val_abs = abs(val)
        z_str = power == 0 ? "" : "z^-$power"
        return " $sign_str $(@sprintf("%.4f", val_abs))$z_str"
    end
    
    # 拼接字符串
    num_str = @sprintf("%.4f", b_sec[1]) * fmt_term(b_sec[2], 1) * fmt_term(b_sec[3], 2)
    den_str = @sprintf("%.4f", a_sec[1]) * fmt_term(a_sec[2], 1) * fmt_term(a_sec[3], 2)
    
    println("子系统 H$i(z):")
    println("        $num_str")
    println("  ----------------------------")
    println("        $den_str")
    println("")
end