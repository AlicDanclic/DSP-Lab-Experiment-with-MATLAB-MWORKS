using TyPlot
using Statistics
using Random

# ==========================================
# 1) 5 点中值滤波器（边界：缩窗）
# ==========================================
"""
    median_filter(signal, L)

一维中值滤波器（窗口长度 L 必须为奇数）。
边界处采用“缩窗”（不补零、不镜像），避免边缘失真。
"""
function median_filter(signal::AbstractVector{T}, L::Int) where {T<:Real}
    @assert L > 0 && isodd(L) "窗口长度 L 必须为正奇数（例如 5）"

    N = length(signal)
    y = Vector{Float64}(undef, N)
    half = L ÷ 2

    @inbounds for i in 1:N
        a = max(1, i - half)
        b = min(N, i + half)
        @views y[i] = median(signal[a:b])
    end
    return y
end

# ==========================================
# 2) 题目数据：s[n] 与噪声 d[n]
# ==========================================
N = 40
n = collect(0:N-1)                 # 用 collect 避免绘图函数对 Range 支持不一致

# s[n] = 3 * n * (0.8)^n
s = 3.0 .* n .* (0.8) .^ n

# d[n]：幅度 0.6 的随机序列（均匀分布在 [-0.6, 0.6]）
noise_amp = 0.6
Random.seed!(2026)                 # 固定随机种子：每次运行结果一致，便于验收/调参
d = noise_amp .* (2 .* rand(N) .- 1)

# 受干扰序列 x[n]
x = s .+ d

# ==========================================
# 3) 5 点中值滤波
# ==========================================
L = 5
y = median_filter(x, L)

# ==========================================
# 4) 绘图：按题意分别画 40 点“受干扰序列”和“滤波输出”
#    （并额外用虚线叠加 s[n] 方便对比，噪声会更明显）
# ==========================================
figure("Median Filter (L=5)", figsize=(10, 8))

# (1) 受干扰序列
subplot(4, 1, 1)
# TyPlot 的 basefmt 不支持 " "（空格式）。这里先画 stem，然后把 baseline 隐藏掉。
h1 = stem(n, x, markerfmt="bo", linefmt="b-", label="x[n]=s[n]+d[n]")
try
    c1 = (h1 isa AbstractVector && length(h1) == 1) ? h1[1] : h1
    if PyCall.pyhasattr(c1, "baseline")
        c1.baseline.set_visible(false)
    end
catch
end
plot(n, s, "g--", linewidth=2, alpha=0.8, label="s[n] (无噪声)")

title("输入信号")
ylabel("幅度")
legend(loc="best")
grid(true)


# (1) 受干扰序列
subplot(4, 1,2)
# TyPlot 的 basefmt 不支持 " "（空格式）。这里先画 stem，然后把 baseline 隐藏掉。
h2 = stem(n, x, markerfmt="bo", linefmt="b-", label="x[n]=s[n]+d[n]")
try
    c2 = (h2 isa AbstractVector && length(h2) == 1) ? h2[1] : h2
    if PyCall.pyhasattr(c2, "baseline")
        c2.baseline.set_visible(false)
    end
catch
end

# 额外把噪声画出来（看不清噪声时非常有用）
plot(n, d, "k:", linewidth=1.5, alpha=0.8, label="d[n] (噪声)")

title("长度为 40 的受干扰序列（噪声幅度 ±0.6）")
ylabel("幅度")
legend(loc="best")
grid(true)


# (3) 中值滤波输出
subplot(4, 1, 3)
h3 = stem(n, y, markerfmt="ro", linefmt="r-", label="信号+噪声")
try
    c3 = (h3 isa AbstractVector && length(h3) == 1) ? h3[1] : h3
    if PyCall.pyhasattr(c3, "baseline")
        c3.baseline.set_visible(false)
    end
catch
end
plot(n, s+d, "k:", linewidth=2, alpha=0.8, label="信号+噪声")

title("信号+噪声")
xlabel("样本序号 n")
ylabel("幅度")
legend(loc="best")
grid(true)

# (2) 中值滤波输出
subplot(4, 1, 4)
h4 = stem(n, y, markerfmt="ro", linefmt="r-", label="y[n] (5点中值滤波输出)")
try
    c4 = (h4 isa AbstractVector && length(h4) == 1) ? h4[1] : h4
    if PyCall.pyhasattr(c4, "baseline")
        c4.baseline.set_visible(false)
    end
catch
end
plot(n, y, "k:", linewidth=2, alpha=0.8, label="y[n] (5点中值滤波输出)")

title("5 点中值滤波器输出")
xlabel("样本序号 n")
ylabel("幅度")
legend(loc="best")
grid(true)

try
    gcf().tight_layout()
catch
end


# 如果你只想严格按题意“分别绘制”且不叠加 s[n] 和 d[n]：
# - 把上面两个 subplot 里的 plot(n, s, ...) 和 plot(n, d, ...) 注释掉即可。