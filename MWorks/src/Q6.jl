import Pkg

# 自动检查并加载 TyPlot 绘图库
if Base.find_package("TyPlot") === nothing
    println("正在安装 TyPlot...")
    Pkg.add("TyPlot")
end

using TyPlot
using LinearAlgebra

println("=== 切比雪夫 I 型模拟低通滤波器设计 (MWorks/Julia) ===")

# ==========================================
# 1. 辅助函数：处理用户输入
# ==========================================
function get_user_input(prompt::String, parse_type::Type, validator::Function=(x->true))
    while true
        print(prompt)
        flush(stdout) # 确保提示符立即显示
        try
            str = strip(readline())
            if isempty(str) continue end
            val = parse(parse_type, str)
            if validator(val)
                return val
            else
                println("输入数值超出合理范围，请重新输入。")
            end
        catch
            println("输入格式无效，请输入一个有效的 $(parse_type)。")
        end
    end
end

# ==========================================
# 2. 核心算法：设计滤波器
# ==========================================
"""
    design_chebyshev_lowpass(N, fc_3db, Rp_db)

参数:
  - N: 阶数
  - fc_3db: 用户输入的 3dB 截止频率
  - Rp_db: 通带波纹 (dB)
返回:
  - poles: S平面极点
  - numerator: 分子系数 (增益)
  - fp_edge: 计算出的通带边缘频率 (在此频率处衰减量为 Rp)
"""
function design_chebyshev_lowpass(N::Int, fc_3db::Float64, Rp_db::Float64)
    # 1. 计算波纹因子 epsilon
    epsilon = sqrt(10^(Rp_db / 10.0) - 1.0)
    
    # 2. 计算缩放因子
    # 切比雪夫滤波器的 3dB 频率与通带边缘频率的关系：
    # ω_c = ωp * cosh(1/N * acosh(1/ε))
    scaling_factor = cosh((1.0 / N) * acosh(1.0 / epsilon))
    
    # 用户输入的 fc_3db 就是 -3dB 频率
    wp_3db = 2 * pi * fc_3db  # 3dB 频率的角频率
    wp_edge = wp_3db / scaling_factor  # 通带边缘的角频率
    fp_edge = wp_edge / (2 * pi)  # 通带边缘频率 (Hz)
    
    # 3. 计算 S 平面极点 (分布在椭圆上)
    mu = asinh(1.0 / epsilon) / N
    poles = ComplexF64[]
    
    for k in 1:N
        # 切比雪夫极点角度公式
        theta = (2*k - 1) * pi / (2*N)
        
        sigma = -sinh(mu) * sin(theta)
        omega = cosh(mu) * cos(theta)
        
        # 使用通带边缘频率进行缩放
        s_k = wp_edge * (sigma + im * omega)
        push!(poles, s_k)
    end
    
    # 4. 计算归一化增益 (分子系数)
    denom_at_0 = real(prod(-poles))
    
    if N % 2 == 1  # 奇数阶
        numerator = denom_at_0  # 直流增益 = 1
    else  # 偶数阶
        numerator = denom_at_0 / sqrt(1 + epsilon^2)
    end
    
    return poles, numerator, fp_edge
end

# 计算复数频率响应
function calc_complex_response(poles, numerator, f_array)
    w_array = 2 * pi .* f_array
    s_array = im .* w_array
    # 传递函数 H(s) = K / product(s - pk)
    H = map(s -> numerator / prod(s .- poles), s_array)
    return H
end

# ==========================================
# 3. 主程序
# ==========================================
function main()
    try
        # --- 获取参数 ---
        N = get_user_input("请输入滤波器阶数 N (整数, >0): ", Int, x -> x > 0)
        fc_input = get_user_input("请输入 3-dB 截止频率 fc (Hz, >0): ", Float64, x -> x > 0)
        Rp = get_user_input("请输入通带波纹 Rp (dB, >0): ", Float64, x -> x > 0)

        println("\n-------------------------------------------")
        println("正在计算滤波器参数...")
        println("  阶数 (N)       : $N")
        println("  3dB 截止频率   : $fc_input Hz")
        println("  通带波纹       : $Rp dB")
        println("-------------------------------------------")

        # --- 设计滤波器 ---
        poles, numerator, fp_calc = design_chebyshev_lowpass(N, fc_input, Rp)
        
        println("计算完成:")
        println("  通带边缘频率 (fp): $(round(fp_calc, digits=2)) Hz (波纹终止点)")
        
        # --- 准备绘图数据 ---
        # 1. 全景数据 (0 到 3倍截止频率)
        f_full = range(0, stop=3*fc_input, length=1000)
        H_full = calc_complex_response(poles, numerator, f_full)
        mag_full_db = 20 .* log10.(abs.(H_full))

        # 2. 细节数据 (0 到 1.2倍通带边缘，用于观察波纹)
        f_zoom = range(0, stop=1.2*fp_calc, length=1000)
        H_zoom = calc_complex_response(poles, numerator, f_zoom)
        mag_zoom_db = 20 .* log10.(abs.(H_zoom))

        # --- 绘图 ---
        println("\n正在绘制图形...")
        
        # 子图 1: 全景响应
        subplot(2, 1, 1)
        plot(f_full, mag_full_db, "b-", linewidth=2, label="增益响应")
        hold("on")
        # 标记 3dB 点
        plot([fc_input, fc_input], [-100, 10], "r--", linewidth=1.5, label="3dB 截止频率")
        plot([0, maximum(f_full)], [-Rp, -Rp], "g:", linewidth=1.5, label="波纹界限 (-$(Rp)dB)")
        
        title("切比雪夫 I 型滤波器: 整体幅频响应")
        ylabel("增益 (dB)")
        grid("on")
        legend()
        
        # 计算最小增益值，但确保不会低于 -60dB
        min_gain = minimum(mag_full_db)
        y_lower = max(-60.0, min_gain)
        y_upper = Rp + 2
        ylim(y_lower, y_upper) # 限制Y轴范围以便观察衰减

        # 子图 2: 通带细节放大
        subplot(2, 1, 2)
        plot(f_zoom, mag_zoom_db, "b-", linewidth=2, label="通带波纹")
        hold("on")
        # 绘制 0dB 和 -Rp dB 参考线
        plot([0, maximum(f_zoom)], [0, 0], "k--", linewidth=1)
        plot([0, maximum(f_zoom)], [-Rp, -Rp], "g--", linewidth=1, label="波纹下限")
        
        title("细节放大: 通带内的波纹特性")
        xlabel("频率 (Hz)")
        ylabel("增益 (dB)")
        grid("on")
        
        # 根据滤波器阶数设置合适的Y轴范围
        if N % 2 == 1
            # 奇数阶：DC处增益为0dB
            ylim(-Rp*1.2, 0.5)
        else
            # 偶数阶：DC处增益为-Rp dB
            ylim(-Rp*1.5, 0.5)
        end

        println("绘图完成。")
        println("请查看绘图窗口：")
        println("  - 上图展示了整体的滤波衰减效果。")
        println("  - 下图放大了通带部分，您可以清晰地看到波纹抖动。")

    catch e
        println("\n发生错误: $e")
        println("错误类型: $(typeof(e))")
        println("请检查参数设置是否正确。")
    end
end

# 运行主程序
main()