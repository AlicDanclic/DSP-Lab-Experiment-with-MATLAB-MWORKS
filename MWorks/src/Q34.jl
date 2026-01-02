using LinearAlgebra
using TyPlot

# coeffs = [a0, a1, ..., an] 表示 a0 + a1*x + ... + an*x^n
function poly_roots(coeffs::AbstractVector{<:Real})
    n = length(coeffs) - 1
    @assert n >= 1 "多项式阶数必须≥1"
    an = coeffs[end]
    @assert abs(an) > 1e-12 "最高次系数不能为0"

    c = reverse(coeffs) ./ an  # [1, a_{n-1}/an, ..., a0/an]

    C = zeros(ComplexF64, n, n)
    C[1, :] .= -c[2:end]
    for i in 2:n
        C[i, i-1] = 1.0
    end
    eigvals(C)
end

# ---------------------------
# 题目：H(z)=B(z^-1)/A(z^-1)
# 乘 z^4 得到 z 平面多项式：
# Bz(z)=z^4 -0.2 z^3 +0.5 z^2 +2 z -0.6
# Az(z)=z^4 +3.2 z^3 +1.5 z^2 -0.8 z +1.4
# ---------------------------

# 注意：poly_roots 需要的是“升幂”系数 [常数项, z^1, z^2, ... , z^4]
Bz = [-0.6,  2.0, 0.5, -0.2, 1.0]
Az = [ 1.4, -0.8, 1.5,  3.2, 1.0]

zzeros = poly_roots(Bz)
zpoles = poly_roots(Az)

println("Zeros (z-plane) = ", round.(zzeros, digits=4))
println("Poles (z-plane) = ", round.(zpoles, digits=4))

# ===== 绘图：不使用 hold，一次 plot 叠加三组数据 =====
theta = range(0, 2π, length=400)
ucx, ucy = cos.(theta), sin.(theta)

# 为了让零点/极点“只有标记不连线”，用 NaN 打断折线
function nansep(v::AbstractVector{<:Real})
    out = Vector{Float64}(undef, 2length(v))
    out[1:2:end] .= v
    out[2:2:end] .= NaN
    out
end

zx = nansep(real.(zzeros)); zy = nansep(imag.(zzeros))
px = nansep(real.(zpoles)); py = nansep(imag.(zpoles))

figure("Question 34: Pole-Zero Plot")

# 一次画三组：单位圆(k--)、零点(bo)、极点(rx)
plot(ucx, ucy, "k--",
     zx,  zy,  "bo",
     px,  py,  "rx",
     linewidth=1, markersize=9)

title("Pole-Zero Plot of the System")
xlabel("Real Axis")
ylabel("Imaginary Axis")
grid("on")
axis("equal")

# 范围
vals = vcat(real.(zzeros), imag.(zzeros), real.(zpoles), imag.(zpoles), [-1.0, 1.0])
m = maximum(abs.(vals)) + 0.5
xlim([-m, m]); ylim([-m, m])

# 图例（如果这一行不兼容，就删掉）
legend(["Unit Circle", "Zeros", "Poles"])