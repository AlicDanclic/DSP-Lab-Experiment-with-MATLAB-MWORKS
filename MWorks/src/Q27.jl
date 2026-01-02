using LinearAlgebra
using Printf

# =========================================================
# 1) 求多项式根（降幂系数）：c0*z^n + c1*z^(n-1) + ... + cn
# =========================================================
function roots_desc(coeffs::AbstractVector{<:Real})
    n = length(coeffs) - 1
    @assert n >= 1 "多项式阶数必须≥1"
    c0 = coeffs[1]
    @assert abs(c0) > 1e-12 "最高次系数不能为0"

    c = coeffs ./ c0  # 归一化，使最高次系数=1：z^n + c1 z^(n-1)+...+cn
    C = zeros(ComplexF64, n, n)

    # companion matrix:
    # [ 0 0 ... 0 -cn
    #   1 0 ... 0 -c(n-1)
    #   0 1 ... 0 -c(n-2)
    #   ...
    #   0 0 ... 1 -c1 ]
    for i in 2:n
        C[i, i-1] = 1.0
    end
    C[:, n] .= -reverse(c[2:end])  # [-cn, ..., -c1]

    return eigvals(C)
end

# =========================================================
# 2) 共轭配对：输出每个section的根(1个或2个)
#    - 复根按 conjugate 配对（保证二阶节实系数）
#    - 实根两两配对，若剩一个则形成一阶节
# =========================================================
function pair_conjugates(r::Vector{ComplexF64}; tol=1e-6)
    r = copy(r)
    used = falses(length(r))
    pairs = Vector{Vector{ComplexF64}}()

    # 先处理复根（找共轭）
    for i in eachindex(r)
        used[i] && continue
        if abs(imag(r[i])) > tol
            j = findfirst(j -> !used[j] && j != i && abs(r[j] - conj(r[i])) < 1e-4, eachindex(r))
            j === nothing && error("未找到共轭根：$(r[i])（可能 tol 太小或根计算异常）")
            push!(pairs, [r[i], r[j]])
            used[i] = true
            used[j] = true
        end
    end

    # 再处理实根
    real_roots = ComplexF64[]
    for i in eachindex(r)
        used[i] && continue
        if abs(imag(r[i])) <= tol
            push!(real_roots, ComplexF64(real(r[i]), 0.0))
            used[i] = true
        end
    end
    sort!(real_roots, by=x->real(x))

    i = 1
    while i <= length(real_roots)
        if i == length(real_roots)
            push!(pairs, [real_roots[i]])        # 一阶节
            i += 1
        else
            push!(pairs, [real_roots[i], real_roots[i+1]])  # 二阶节
            i += 2
        end
    end

    return pairs
end

# =========================================================
# 3) 由 z 平面根构造 z^-1(=q) 形式的节系数
#    若节含根 z1,z2：
#       (1 - z1 q)(1 - z2 q) = 1 -(z1+z2)q + (z1 z2) q^2
#    若节含单根 z1：
#       (1 - z1 q)
# =========================================================
function section_coeff_from_zroots(zroots::Vector{ComplexF64})
    if length(zroots) == 2
        z1, z2 = zroots[1], zroots[2]
        b0 = 1.0
        b1 = -real(z1 + z2)
        b2 =  real(z1 * z2)
        return [b0, b1, b2]  # 升幂：b0 + b1 q + b2 q^2
    elseif length(zroots) == 1
        z1 = zroots[1]
        b0 = 1.0
        b1 = -real(z1)
        return [b0, b1]      # 一阶：b0 + b1 q
    else
        return [1.0]         # 空节（不会用到）
    end
end

# =========================================================
# 4) 多项式卷积（升幂系数：常数项在前）
# =========================================================
function conv(a::Vector{Float64}, b::Vector{Float64})
    y = zeros(Float64, length(a) + length(b) - 1)
    for i in eachindex(a)
        for j in eachindex(b)
            y[i+j-1] += a[i] * b[j]
        end
    end
    return y
end

# =========================================================
# 5) 主程序：题目系统
#    H(z)= (2+5z^-1+z^-2-3z^-3+4z^-4+6z^-5) / (1+3z^-1-5z^-2+2z^-3-4z^-4+3z^-5)
# =========================================================
b = [2.0, 5.0, 1.0, -3.0, 4.0, 6.0]   # B(q), q=z^-1, 升幂
a = [1.0, 3.0, -5.0, 2.0, -4.0, 3.0]  # A(q), 升幂

# z 平面零极点来自：z^5 B(1/z) 与 z^5 A(1/z)
# 它们的降幂系数正好是 [b0,b1,...,b5] 与 [a0,a1,...,a5]（对应 z^5 + ... + 常数）
Bz_desc = b               # 2 z^5 +5 z^4 +... +6
Az_desc = a               # 1 z^5 +3 z^4 +... +3

zzeros = roots_desc(Bz_desc)
zpoles = roots_desc(Az_desc)

println("z-plane zeros:")
println.(round.(zzeros, digits=6))
println("\nz-plane poles:")
println.(round.(zpoles, digits=6))

# 配对成级联节
zero_pairs = pair_conjugates(zzeros)
pole_pairs = pair_conjugates(zpoles)

# 级联节数量（对齐：不足的补1）
S = max(length(zero_pairs), length(pole_pairs))
while length(zero_pairs) < S
    push!(zero_pairs, ComplexF64[])  # 表示“无零点”的节 => 分子=1
end
while length(pole_pairs) < S
    push!(pole_pairs, ComplexF64[])  # 表示“无极点”的节 => 分母=1
end

# 总增益（常数项比值）
k = b[1] / a[1]

println("\n==============================")
println("级联型结构：H(z) = k * Π H_i(z)")
@printf("k = %.6f\n", k)
println("==============================\n")

# 构造并打印每个子系统
num_sections = Vector{Vector{Float64}}()
den_sections = Vector{Vector{Float64}}()

for i in 1:S
    bn = section_coeff_from_zroots(zero_pairs[i])
    ad = section_coeff_from_zroots(pole_pairs[i])

    push!(num_sections, bn)
    push!(den_sections, ad)

    # 打印表达式（按 z^-1）
    if length(bn) == 3
        @printf("H_%d(z) 分子:  1 %+ .6f z^-1 %+ .6f z^-2\n", i, bn[2], bn[3])
    elseif length(bn) == 2
        @printf("H_%d(z) 分子:  1 %+ .6f z^-1\n", i, bn[2])
    else
        @printf("H_%d(z) 分子:  1\n", i)
    end

    if length(ad) == 3
        @printf("      分母:  1 %+ .6f z^-1 %+ .6f z^-2\n\n", ad[2], ad[3])
    elseif length(ad) == 2
        @printf("      分母:  1 %+ .6f z^-1\n\n", ad[2])
    else
        @printf("      分母:  1\n\n")
    end
end

function check_reconstruct(num_sections, den_sections, k, b, a)
    # 6) 重构校验：把各节卷起来，检查是否回到原 b,a
    b_rec = [1.0]
    a_rec = [1.0]
    for i in 1:length(num_sections)
        b_rec = conv(b_rec, num_sections[i])
        a_rec = conv(a_rec, den_sections[i])
    end
    b_rec .*= k

    while length(b_rec) < length(b)
        push!(b_rec, 0.0)
    end
    while length(a_rec) < length(a)
        push!(a_rec, 0.0)
    end

    println("===== 重构校验（应与原系数一致）=====")
    println("原 b = ", b)
    println("重构b = ", round.(b_rec, digits=6))
    println("原 a = ", a)
    println("重构a = ", round.(a_rec, digits=6))

    err_b = maximum(abs.(b_rec[1:length(b)] .- b))
    err_a = maximum(abs.(a_rec[1:length(a)] .- a))
    @printf("max|Δb| = %.3e,  max|Δa| = %.3e\n", err_b, err_a)
end

check_reconstruct(num_sections, den_sections, k, b, a)