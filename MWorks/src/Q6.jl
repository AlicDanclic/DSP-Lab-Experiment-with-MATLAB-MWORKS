import Pkg
# 检查并添加 TyPlot
if Base.find_package("TyPlot") === nothing
    println("正在安装 TyPlot ...")
    Pkg.add("TyPlot")
end

using TyPlot
using LinearAlgebra

println("=== 切比雪夫 I 型模拟低通滤波器设计 ===")

# ==========================================
# 1. 键盘输入参数
# ==========================================
function get_valid_input(prompt::String, parse_func::Function)
    while true
        print(prompt)
        flush(stdout)
        try
            input_str = strip(readline())
            if isempty(input_str) continue end
            return parse_func(input_str)
        catch
            println("输入格式不正确，请重新输入。")
        end
    end
end

try
    # 获取参数
    global N = get_valid_input("请输入滤波器阶数 N (整数, 例如 4): ", x -> parse(Int, x))
    global fc = get_valid_input("请输入截止频率 fc (Hz, 例如 1000): ", x -> parse(Float64, x))
    global Rp = get_valid_input("请输入通带波纹 Rp (dB, 例如 0.5): ", x -> parse(Float64, x))

    println("\n正在设计: N = $N, fc = $fc Hz, Rp = $Rp dB")

catch e
    println("\n错误: $e")
    exit()
end

# ==========================================
# 2. 计算极点 (S平面)
# ==========================================
# 切比雪夫 I 型滤波器的极点位于 S 平面的椭圆上
# 这里的 fc 被视为通带边缘频率 (Passband Edge Frequency)

wc = 2 * pi * fc
epsilon = sqrt(10^(Rp/10) - 1) # 波纹因子
mu = asinh(1/epsilon) / N      # 辅助参数

poles = ComplexF64[]
for k in 1:N
    # 角度定义
    # theta_k = pi/2 + (2k-1)/(2N) * pi
    # 注意: 有些定义 k 从 0 到 N-1，这里 k 从 1 到 N，公式对应调整
    theta = pi/2 + ((2*k - 1) * pi) / (2 * N)
    
    # 极点公式: s_k = wc * (-sinh(mu)sin(theta) + j cosh(mu)cos(theta))
    real_part = -sinh(mu) * sin(theta)
    imag_part = cosh(mu) * cos(theta)
    s_k = wc * (real_part + im * imag_part)
    
    push!(poles, s_k)
end

println("\n计算得到的极点 (S平面):")
for (i, p) in enumerate(poles)
    println("p$i: $(round(real(p), digits=2)) + $(round(imag(p), digits=2))j")
end

# ==========================================
# 3. 计算频率响应
# ==========================================
f_plot = range(0, stop=3*fc, length=1000)
w_plot = 2 * pi .* f_plot

# 确定分子系数 K
# 切比雪夫 I 型滤波器在 w=0 处的幅值为:
# - 如果 N 是奇数: |H(0)| = 1 (0 dB)
# - 如果 N 是偶数: |H(0)| = 1 / sqrt(1 + epsilon^2) (-Rp dB)

# 计算分母多项式在 s=0 处的值 (即所有极点的积的模)
denom_at_0 = abs(prod(-poles)) # s=0 -> product(-pk)

target_dc_gain = (N % 2 == 1) ? 1.0 : (1.0 / sqrt(1 + epsilon^2))
numerator = target_dc_gain * denom_at_0

gains = Float64[]

for w in w_plot
    s = im * w
    # H(s) = K / product(s - pk)
    denominator = prod(s .- poles)
    H = numerator / denominator
    push!(gains, abs(H))
end

gains_db = 20 .* log10.(gains)

# ==========================================
# 4. 绘图
# ==========================================
println("\n正在绘制增益响应...")

plot(f_plot, gains_db, "b-", linewidth=2, label="Gain Response")
hold("on")

# 绘制通带波纹下限
plot([0, fc], [-Rp, -Rp], "g--", label="Passband Ripple (-$Rp dB)")

# 标记截止频率点
plot([fc, fc], [-60, 5], "r:", label="Cutoff Frequency")

title("Chebyshev Type I Analog Filter (N=$N, fc=$(Int(fc))Hz, Rp=$(Rp)dB)")
xlabel("Frequency (Hz)")
ylabel("Magnitude (dB)")
grid("on")
legend()

# 调整 Y 轴范围以便观察波纹
ylim(min(-40, -Rp-10), max(5, Rp+2))

println("绘图完成。")