using LinearAlgebra

# ========= 工具：求根（输入为“降幂”系数：c0*x^n + c1*x^(n-1)+...+cn）
function roots_desc(coeffs::AbstractVector{<:Real})
    n = length(coeffs) - 1
    @assert n >= 1
    lead = coeffs[1]
    @assert abs(lead) > 1e-12
    c = coeffs ./ lead               # 归一化成首项为1

    C = zeros(ComplexF64, n, n)
    for i in 1:n-1
        C[i+1, i] = 1.0
    end
    C[:, n] .= -c[2:end]
    eigvals(C)
end

# ========= 多项式求值（升幂：a0 + a1*q + ...）
polyval_asc(a, q) = sum(a[i] * q^(i-1) for i in 1:length(a))

# ========= 多项式导数（升幂系数 -> 升幂系数）
function polyder_asc(a::AbstractVector{<:Real})
    n = length(a) - 1
    n == 0 && return [0.0]
    d = zeros(Float64, n)
    for i in 2:length(a)
        d[i-1] = (i-1) * a[i]
    end
    d
end

# ========= MATLAB residuez 等价实现：b,a 为 z^-1(=q) 的升幂系数 [b0,b1,...]
# 返回 r,p,k，使 H(z)=k + Σ r_i/(1 - p_i z^-1)
function residuez_like_matlab(b::Vector{Float64}, a::Vector{Float64})
    nb = length(b)-1
    na = length(a)-1
    @assert na >= 1

    # 这里只写够用版本：若 nb==na，则 k 为常数（最高次项比值）
    # 你的题刚好 nb==na==5
    @assert nb == na "如需处理 nb≠na，可再扩展做多项式长除法"
    k = b[end] / a[end]                  # 注意：最高次项比值（不是 b[1]/a[1]）
    b_rem = b .- k .* a                  # 余式（仍是升幂系数）

    # q=z^-1 域极点：A(q)=0 的根
    q_poles = roots_desc(reverse(a))     # reverse(a) 变成 q 的降幂系数
    p = 1.0 ./ q_poles                   # z 域极点

    # 留数：res_q = B_rem(q_i)/A'(q_i)，再换成 MATLAB 形式 r_i = -p_i * res_q
    a_der = polyder_asc(a)
    r = similar(p)
    for i in eachindex(q_poles)
        qi = q_poles[i]
        res_q = polyval_asc(b_rem, qi) / polyval_asc(a_der, qi)
        r[i] = -(p[i]) * res_q
    end

    # MATLAB 的 k 会把由换元产生的常数项也吸收进去；
    # 对于 r/(1-p z^-1)=r + (r p)/(z-p)，常数项是 r，所以把 Σr 加到 k 里更接近 MATLAB 输出
    k = k + sum(real.(r))  # 系统系数全实，最终 k 应为实数
    return r, p, k
end

# ========= 题目系数
b = [2.0, 5.0, 1.0, -3.0, 4.0, 6.0]
a = [1.0, 3.0, -5.0, 2.0, -4.0, 3.0]

r, p, k = residuez_like_matlab(b, a)

println("Parallel Form I (complex):")
println("k = ", k)
for i in eachindex(p)
    println("r[$i] = ", r[i], " , p[$i] = ", p[i])
end

println("\nParallel Form II (real SOS):")
tol = 1e-6
used = falses(length(p))
for i in eachindex(p)
    used[i] && continue
    if abs(imag(p[i])) < tol
        println("  section: ", real(r[i]), " / (1 - ", real(p[i]), " z^-1)")
        used[i] = true
    else
        # 找共轭
        j = findfirst(j -> !used[j] && abs(p[j] - conj(p[i])) < 1e-5, eachindex(p))
        @assert j !== nothing "未找到共轭极点，请检查 tol"
        # 合并为二阶节（实系数）
        b0 = 2*real(r[i])
        b1 = -2*real(r[i]*conj(p[i]))
        a1 = -2*real(p[i])
        a2 = abs(p[i])^2
        println("  section: (", b0, " + ", b1, " z^-1) / (1 + ", a1, " z^-1 + ", a2, " z^-2)")
        used[i] = true
        used[j] = true
    end
end