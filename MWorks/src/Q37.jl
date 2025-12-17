using TyPlot        # Syslab 原生绘图
using LinearAlgebra
using Printf

# ==========================================
#  1. 巴特沃兹滤波器设计算法
# ==========================================

function design_butterworth_analog(N, Wc)
    # 计算极点位置
    # 极点分布在半径为 Wc 的左半圆上
    # 角度公式: theta_k = pi * (2k + N - 1) / (2N), k = 1...N
    
    poles = ComplexF64[]
    for k in 1:N
        theta = pi * (2*k + N - 1) / (2*N)
        s = Wc * exp(im * theta)
        push!(poles, s)
    end
    
    # 根据极点构建分母多项式系数 D(s) = (s-p1)(s-p2)...
    # 我们使用多项式乘法展开
    # 初始化多项式为 [1.0] (对应常数 1)
    den_poly = [1.0 + 0.0im] 
    
    for p in poles
        # 多项式乘法: (coeff...)*s - (coeff...)*p
        # 相当于卷积: conv(den_poly, [1, -p])
        new_poly = zeros(ComplexF64, length(den_poly) + 1)
        # 移位加
        new_poly[1:end-1] += den_poly    # 对应 * s
        new_poly[2:end]   -= den_poly * p # 对应 * -p
        den_poly = new_poly
    end
    
    # 理论上巴特沃兹多项式系数全是实数，取实部即可
    return real.(den_poly)
end

# ==========================================
#  2. 主程序
# ==========================================

N_order = 4
Wc = 1.0

println("正在设计 4 阶巴特沃兹滤波器 (Wc=1)...")
den_coeffs = design_butterworth_analog(N_order, Wc)

# 打印传递函数
println("\n========================================")
println("       设计结果: 传递函数 H(s)")
println("========================================")
println("H(s) = 1 / D(s)")
println("D(s) 系数 (从高阶 s^4 到 s^0):")
println(round.(den_coeffs, digits=4))
println("\n数学表达式:")
print("H(s) = 1 / ( s^4")
@printf(" + %.4fs^3", den_coeffs[2])
@printf(" + %.4fs^2", den_coeffs[3])
@printf(" + %.4fs",   den_coeffs[4])
@printf(" + %.4f )",  den_coeffs[5])
println("\n")

# ==========================================
#  3. 计算增益响应并绘图
# ==========================================

# 定义频率范围 (模拟频率 rad/s)
w_axis = range(0, 3.0, length=500)

# 直接利用模平方公式计算幅度，这样最精确
# |H(jw)| = 1 / sqrt(1 + (w/Wc)^(2N))
mag_response = [1.0 / sqrt(1 + (w/Wc)^(2*N_order)) for w in w_axis]

# 转换为 dB
gain_db = 20 * log10.(mag_response)

figure("Question 37: Butterworth Response")

# 1. 绘制增益曲线
plot(w_axis, gain_db, "b-", linewidth=2, label="Gain Response")
hold("on")

# 2. 标记 3-dB 截止点
# 在 w=1 处，增益应为 -3.01 dB
plot([1.0], [-3.0103], "ro", markersize=8, label="Cutoff (1 rad/s, -3dB)")

# 3. 绘制辅助线
plot([0, 1.0], [-3.01, -3.01], "r--", linewidth=1) # 横向虚线
plot([1.0, 1.0], [-50, -3.01], "r--", linewidth=1) # 纵向虚线

title("Gain Response of 4th-Order Butterworth Low-Pass Filter")
xlabel("Frequency (rad/s)")
ylabel("Gain (dB)")
grid("on")
legend("on")

# 设置显示范围
ylim([-40, 5])
xlim([0, 3])

hold("off")