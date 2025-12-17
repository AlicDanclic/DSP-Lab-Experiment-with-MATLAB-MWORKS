using LinearAlgebra
using TyPlot # 使用 Syslab 原生绘图

# =========================================================================
#  第一部分：FIR 滤波器设计核心算法 (窗函数法)
# =========================================================================

function design_multiband_fir()
    # 1. 确定参数
    # 过渡带宽度 delta_w = 0.03pi
    # 海明窗估算阶数: N ≈ 6.6pi / delta_w ≈ 6.6 / 0.03 = 220
    # 我们选择长度 L = 221 (偶对称，线性相位)
    N_order = 220
    L = N_order + 1
    tau = N_order / 2  # 中心点
    
    # 截止频率 (取过渡带中心)
    wc1 = 0.335 * pi  # (0.32 + 0.35)/2
    wc2 = 0.665 * pi  # (0.65 + 0.68)/2
    
    # 2. 计算理想脉冲响应 h_d[n]
    # 理想幅频特性可以看作三个矩形频带的叠加，或者全通+低通的组合。
    # 这里我们利用线性性质分解：
    # FullBand(0.8) - LowPass_wc2(0.6) + LowPass_wc1(0.4)
    # 解释:
    #   0 ~ wc1: 0.8 - 0.6 + 0.4 = 0.6 (满足)
    #   wc1 ~ wc2: 0.8 - 0.6 + 0 = 0.2 (满足)
    #   wc2 ~ pi:  0.8 - 0 + 0 = 0.8 (满足)
    
    h = zeros(Float64, L)
    for i in 0:L-1
        n_shift = i - tau
        
        # 处理 n=0 (L'Hopital 极限)
        if abs(n_shift) < 1e-9
            # h[0] = 0.8 - 0.6*(wc2/pi) + 0.4*(wc1/pi)
            val = 0.8 - 0.6*(wc2/pi) + 0.4*(wc1/pi)
        else
            # h[n] = 0.8*delta[n] - 0.6*sin(wc2*n)/(pi*n) + 0.4*sin(wc1*n)/(pi*n)
            # 注意：全通项 0.8*delta[n] 只在 n=0 时存在，这里 n!=0，所以只有 sinc 项
            val = -0.6 * sin(wc2 * n_shift) / (pi * n_shift) + 
                   0.4 * sin(wc1 * n_shift) / (pi * n_shift)
        end
        h[i+1] = val
    end
    
    # 3. 生成海明窗 (Hamming Window)
    # w[n] = 0.54 - 0.46 * cos(2pi*n / (L-1))
    w_win = [0.54 - 0.46 * cos(2 * pi * i / (L - 1)) for i in 0:L-1]
    
    # 4. 加窗得到最终系数
    h_final = h .* w_win
    
    return h_final
end

# 手动计算幅频响应 (DTFT)
function calc_freq_response(h, steps=1000)
    w_vals = range(0, pi, length=steps)
    mag_vals = Float64[]
    
    for w in w_vals
        # H(w) = sum(h[n] * exp(-j*w*n))
        H = 0.0 + 0.0im
        for i in 1:length(h)
            H += h[i] * exp(-im * w * (i-1))
        end
        push!(mag_vals, abs(H))
    end
    return w_vals, mag_vals
end

# =========================================================================
#  第二部分：主程序与绘图
# =========================================================================

println("正在设计 FIR 滤波器 (N=220)...")
h_coeff = design_multiband_fir()

println("正在计算幅频特性...")
w_axis, mag_resp = calc_freq_response(h_coeff)

println("正在绘图...")
figure("Question 32: Multiband FIR Filter")

# 绘制幅频特性
plot(w_axis/pi, mag_resp, "b-", linewidth=1.5, label="Designed Filter")

# 绘制理想的规约线（红色虚线框），方便验证
hold("on")
# 0 ~ 0.32pi: 0.6
plot([0, 0.32], [0.6, 0.6], "r--", linewidth=2, label="Ideal Spec")
# 0.35pi ~ 0.65pi: 0.2
plot([0.35, 0.65], [0.2, 0.2], "r--", linewidth=2)
# 0.68pi ~ 1.0pi: 0.8
plot([0.68, 1.0], [0.8, 0.8], "r--", linewidth=2)

title("Frequency Response of Multiband FIR Filter (Hamming Window)")
xlabel("Normalized Frequency (x pi)")
ylabel("Magnitude")
legend("on")
grid("on")
ylim([0, 1.0]) # 设置纵坐标范围方便观察

hold("off")