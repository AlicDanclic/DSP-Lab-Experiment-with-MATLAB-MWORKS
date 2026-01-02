using LinearAlgebra
using Printf

println("=== Q3: 部分分式展开（对应 MATLAB residuez）===")

# ------------------------------------------------------------
# 1) 题目系数（按 z^-1 升幂：b0 + b1 z^-1 + ...）
# ------------------------------------------------------------
b = [1.0, -0.2, 0.5]                 # 1 - 0.2 z^-1 + 0.5 z^-2
a = [1.0,  3.2, 1.5, -0.8, 1.4]      # 1 + 3.2 z^-1 + 1.5 z^-2 - 0.8 z^-3 + 1.4 z^-4

# ------------------------------------------------------------
# 2) 工具：降幂多项式求根  c0 z^n + c1 z^(n-1)+...+cn
# ------------------------------------------------------------
function roots_desc(c::AbstractVector{<:Real})
    n = length(c) - 1
    @assert n >= 1 "多项式阶数必须≥1"
    c0 = c[1]
    @assert abs(c0) > 1e-12 "最高次系数不能为0"

    cnorm = c ./ c0  # 归一化最高次为1

    C = zeros(ComplexF64, n, n)
    for i in 2:n
        C[i, i-1] = 1.0
    end
    # companion 最后一列 = -[cn, c(n-1), ..., c1]^T
    C[:, n] .= -reverse(cnorm[2:end])
    eigvals(C)
end

# 多项式求值：降幂系数
polyval_desc(c::AbstractVector, z) = begin
    y = 0.0 + 0.0im
    for i in eachindex(c)
        y = y*z + c[i]
    end
    y
end

# 降幂导数系数
function polyder_desc(c::AbstractVector{<:Real})
    n = length(c) - 1
    d = zeros(Float64, n)
    for i in 1:n
        d[i] = c[i] * (n - (i - 1))
    end
    d
end

# ------------------------------------------------------------
# 3) 把 A(z^-1) 转成 z 平面 D(z)=z^M A(z^-1) 并求极点
#    A(z^-1)=1+a1 z^-1+...+aM z^-M
#    => D(z)=z^M + a1 z^(M-1) + ... + aM
#    其降幂系数刚好是 a 本身
# ------------------------------------------------------------
M = length(a) - 1
D_desc = a                    # z^4 + 3.2 z^3 + 1.5 z^2 - 0.8 z + 1.4
p = roots_desc(D_desc)        # z 平面极点

# ------------------------------------------------------------
# 4) 留数计算（对应 H(z)=Σ r_i / (1 - p_i z^-1)）
#
# 令 N(z)=z^M B(z^-1)
# B(z^-1)=b0+b1 z^-1+...+bK z^-K
# => N(z)=b0 z^M + b1 z^(M-1) + ... + bK z^(M-K) + ... (不足补0)
#
# 对简单极点，有公式：
#   r_i = N(p_i) / ( p_i * D'(p_i) )
# ------------------------------------------------------------
# 构造 N(z) 的降幂系数（长度 M+1）
N_desc = zeros(Float64, M+1)
for k in 0:(length(b)-1)
    N_desc[1+k] = b[1+k]  # b0,b1,b2...
end
# N_desc 现在是 [b0,b1,b2,0,0] => z^4 -0.2 z^3 +0.5 z^2

Dder_desc = polyder_desc(D_desc)

r = ComplexF64[]
for pi in p
    ri = polyval_desc(N_desc, pi) / (pi * polyval_desc(Dder_desc, pi))
    push!(r, ri)
end

# ------------------------------------------------------------
# 5) 打印：复数留数形式（与 residuez 同型）
# ------------------------------------------------------------
println("\n--- 极点 p（z 平面）---")
for (i, val) in enumerate(p)
    @printf("p%d = %.6f %+.6fim\n", i, real(val), imag(val))
end

println("\n--- 留数 r（使 H(z)=Σ r/(1-p z^-1) ）---")
for (i, val) in enumerate(r)
    @printf("r%d = %.6f %+.6fim\n", i, real(val), imag(val))
end

println("\n--- 部分分式展开（复数形式）---")
println("H(z) = ")
for i in eachindex(r)
    @printf("  %+.6f %+.6fim / (1 - (%.6f %+.6fim) z^-1)\n",
            real(r[i]), imag(r[i]), real(p[i]), imag(p[i]))
end

# ------------------------------------------------------------
# 6) 可选：把共轭对合并成“实系数二阶项”（更像手写答案）
#    若 (r,p) 与 (conj(r),conj(p))：
#      r/(1-pq) + r*/(1-p*q) =
#      (2Re(r) - 2Re(r*conj(p)) q) / (1 - 2Re(p) q + |p|^2 q^2),  q=z^-1
# ------------------------------------------------------------
function pair_conj_indices(vals::Vector{ComplexF64}; tol=1e-6)
    used = falses(length(vals))
    pairs = Tuple{Int,Int}[]
    singles = Int[]
    for i in eachindex(vals)
        used[i] && continue
        if abs(imag(vals[i])) < tol
            push!(singles, i); used[i]=true
        else
            j = findfirst(j -> !used[j] && j!=i && abs(vals[j]-conj(vals[i])) < 1e-4, eachindex(vals))
            if j === nothing
                push!(singles, i); used[i]=true
            else
                push!(pairs, (i,j)); used[i]=true; used[j]=true
            end
        end
    end
    return pairs, singles
end

pairs, singles = pair_conj_indices(p)

println("\n--- 合并共轭后的实系数形式 ---")
for (i,j) in pairs
    ri, pi = r[i], p[i]
    b0 = 2*real(ri)
    b1 = -2*real(ri*conj(pi))
    a1 = -2*real(pi)
    a2 = abs(pi)^2
    @printf("  (%.6f %+ .6f z^-1) / (1 %+ .6f z^-1 %+ .6f z^-2)\n", b0, b1, a1, a2)
end
for i in singles
    # 实极点对应一阶项
    @printf("  (%.6f) / (1 - (%.6f) z^-1)\n", real(r[i]), real(p[i]))
end

# ------------------------------------------------------------
# 7) 数值校验：在单位圆上采样对比原系统与展开式
#    （包进函数，避免 Mworks soft scope 报错）
# ------------------------------------------------------------
function validate_pf(b, a, r, p)
    function H_from_ba(b,a,ω)
        q = exp(-1im*ω)  # z^-1
        num = sum(b[k]*q^(k-1) for k in 1:length(b))
        den = sum(a[k]*q^(k-1) for k in 1:length(a))
        num/den
    end

    function H_from_rp(r,p,ω)
        q = exp(-1im*ω)
        s = 0.0 + 0.0im
        for i in eachindex(r)
            s += r[i] / (1 - p[i]*q)
        end
        s
    end

    ws = range(0, 2π, length=2048)
    err = 0.0
    for ω in ws
        err = max(err, abs(H_from_ba(b,a,ω) - H_from_rp(r,p,ω)))
    end
    @printf("\n校验：max |H_orig - H_pf| = %.3e\n", err)
end

validate_pf(b, a, r, p)