using TyPlot       # Syslab 原生绘图
using LinearAlgebra

# ==========================================
#  1. 定义 DFT 计算函数 (防止缺少 FFTW 包)
# ==========================================
function my_dft(x)
    N = length(x)
    X = zeros(ComplexF64, N)
    # DFT 公式: X[k] = sum_{n=0}^{N-1} x[n] * exp(-j*2*pi*k*n/N)
    for k in 0:N-1
        sum_val = 0.0 + 0.0im
        for n in 0:N-1
            # 注意：Julia 数组索引从 1 开始，所以取 x[n+1]
            angle = -2 * pi * k * n / N
            sum_val += x[n+1] * exp(im * angle)
        end
        X[k+1] = sum_val
    end
    return X
end

# ==========================================
#  2. 生成信号与计算
# ==========================================
N = 64
n = 0:N-1

# 生成序列 x[n] = sin(25*pi*n / 64)
x = sin.(25 * pi * n / 64)

println("正在计算 64 点 DFT...")
X_k = my_dft(x)

# 计算幅度谱 |X[k]|
mag_X = abs.(X_k)

# ==========================================
#  3. 绘图
# ==========================================
figure("Question 36: DFT Magnitude")

# 使用 stem 图 (火柴梗图) 展示离散频谱
# 也就是画出每一根谱线的高度
k_axis = 0:N-1
stem(k_axis, mag_X, "b-", filled=true, markersize=4, label="|X[k]|")

# 图表修饰
title("Magnitude of 64-point DFT: |X[k]|")
xlabel("Frequency Index k")
ylabel("Magnitude")
grid("on")

# 标记出峰值位置的辅助线 (k=12.5 处)
# 理论中心在 12.5 和 64-12.5=51.5
hold("on")
plot([12.5, 12.5], [0, maximum(mag_X)], "r--", linewidth=1, label="True Freq (k=12.5)")
plot([51.5, 51.5], [0, maximum(mag_X)], "r--", linewidth=1)
legend("on")

hold("off")

println("计算完成。")
println("观察图形：由于频率 k=12.5 不是整数，")
println("频谱能量并未集中在单根谱线上，而是分散在 k=12 和 k=13 周围。")