using TyPlot
using DSP   # 需要安装: 用于 remez 滤波器设计
using FFTW  # 需要安装: 用于计算频率响应 (fft)

# ==============================================================================
# 核心计算函数 1：Hermann 公式估计
# ==============================================================================
function hermann_estimation(fp, fs, dp, ds, Fs)
    # 1. 归一化频率
    Delta_F = (fs - fp) / Fs

    # 2. 预计算对数值
    log_dp = log10(dp)
    log_ds = log10(ds)

    # 3. 计算 D_infinity 系数
    a_dp = 0.005309 * (log_dp^2) + 0.07114 * log_dp - 0.4761
    g_dp = -0.00266 * (log_dp^2) - 0.5941 * log_dp - 0.4278
    D_inf = a_dp * log_ds + g_dp

    # 4. 计算过渡带修正项
    f_correction = 11.01217 + 0.51244 * (log_dp - log_ds)

    # 5. 计算滤波器长度 N
    N_exact = (D_inf - f_correction * (Delta_F^2)) / Delta_F
    
    # 向上取整
    N = ceil(Int, N_exact)
    
    # 阶数 M
    M = N - 1

    return N, M, N_exact
end

# ==============================================================================
# 核心计算函数 2：设计滤波器并计算响应 (修正版)
# ==============================================================================
function design_and_analyze_filter(N, fp, fs, Fs)
    
    # 1. 使用 Parks-McClellan (Remez) 算法设计滤波器
    # 注意：remez 接受的参数是 阶数(Order) = N-1
    M = N - 1
    
    # 定义频带
    # ERROR FIX: bands 必须是一个扁平的向量 [start1, end1, start2, end2]
    # 对应 desired [val1, val2]
    bands = [0.0, fp, fs, Fs/2]
    desired = [1.0, 0.0]
    
    # 权重: 为了更好地满足阻带波纹要求，通常需要给予阻带更高的权重
    # weight = [1.0/dp, 1.0/ds] 或 [1.0, dp/ds]
    # 这里我们根据题目给定的 dp=0.015, ds=0.021 计算权重比
    dp=0.015
    ds=0.021
    w_pass = 1.0 / dp
    w_stop = 1.0 / ds
    weights = [w_pass, w_stop]

    try
        # 调用 remez
        h = remez(M, bands, desired, weight=weights, Hz=Fs)
        
        # 2. 计算频率响应 (使用 FFT)
        n_fft = 8192 # 增加 FFT 点数使曲线更平滑
        H = fft([h; zeros(n_fft - length(h))])
        
        # 3. 生成对应的频率轴
        freq_axis = range(0, Fs, length=n_fft)
        
        # 4. 截取前一半 (0 ~ Nyquist)
        valid_idx = 1:div(n_fft, 2)
        f_plot = freq_axis[valid_idx]
        mag_plot = abs.(H[valid_idx]) # 幅度
        
        return h, f_plot, mag_plot
    catch e
        println("设计滤波器失败: $e")
        println("请检查是否已安装 DSP 包 (import Pkg; Pkg.add(\"DSP\"))")
        return [], [], []
    end
end

# ==============================================================================
# 绘图函数：绘制规格 + 实际响应
# ==============================================================================
function plot_filter_results(fp, fs, dp, ds, Fs, N, M, f_response, mag_response)
    fmax = Fs / 2
    
    figure("FIR Filter Design Result")
    clf()
    
    hold("on")

    # --- 1. 绘制公差框 (Tolerance Boxes) ---
    # 通带框 (蓝色) - 允许范围 [1-dp, 1+dp]
    plot([0, fp], [1+dp, 1+dp], "b-", linewidth=1)
    plot([0, fp], [1-dp, 1-dp], "b-", linewidth=1)
    plot([fp, fp], [1-dp, 1+dp], "b-", linewidth=1) 
    
    # 阻带框 (洋红色) - 允许范围 [0, ds]
    plot([fs, fmax], [ds, ds], "m-", linewidth=1)
    plot([fs, fs], [0, ds], "m-", linewidth=1)

    # --- 2. 绘制实际滤波器响应 ---
    if !isempty(mag_response)
        plot(f_response, mag_response, "k-", linewidth=1.2, label="实际滤波器响应 (N=$N)")
    end

    # --- 3. 辅助线 ---
    plot([0, fmax], [1, 1], "k:", linewidth=0.5) # 1.0 参考线
    plot([0, fmax], [0, 0], "k-", linewidth=0.5) # 0.0 参考线
    
    # 截止频率竖线
    plot([fp, fp], [-0.05, 1+dp+0.05], "k--", alpha=0.3)
    plot([fs, fs], [-0.05, ds+0.05],   "k--", alpha=0.3)

    # --- 4. 标注与设置 ---
    xlim([0, fmax])
    # 稍微放大Y轴范围以便观察波纹
    ylim([-0.05, 1.1]) 
    
    xlabel("频率 (Hz)")
    ylabel("幅度 |H(f)|")
    title("FIR 低通滤波器设计结果 (Hermann 估计阶数 N=$N)")
    grid("on")
    legend("loc", "best")

    # 关键点文字
    text(fp, -0.08, "fp", ha="center")
    text(fs, -0.08, "fs", ha="center")
    
    info_text = "规格:\n通带波纹: $dp\n阻带波纹: $ds"
    text(fmax*0.75, 0.6, info_text, bbox=Dict("facecolor"=>"white", "alpha"=>0.8))

    hold("off")
end

# ==============================================================================
# 主程序
# ==============================================================================

# 1. 定义指标
Fs_val = 5000.0
fp_val = 1500.0
fs_val = 1800.0
dp_val = 0.015
ds_val = 0.021

println("========================================")
println("       Hermann FIR 滤波器设计")
println("========================================")

# 2. 计算阶数
N_est, M_est, N_ex = hermann_estimation(fp_val, fs_val, dp_val, ds_val, Fs_val)

println("Hermann 估算结果:")
println("  精确长度 N_exact = $(round(N_ex, digits=4))")
println("  取整长度 N       = $N_est")
println("  滤波器阶数 M     = $M_est")
println("----------------------------------------")

# 3. 设计滤波器并获取响应数据
h_coef, f_data, mag_data = design_and_analyze_filter(N_est, fp_val, fs_val, Fs_val)

if isempty(h_coef)
    println("警告: 未能生成滤波器数据，仅绘制规格框。")
else
    println("滤波器设计成功。正在绘制响应图...")
end

# 4. 绘图
plot_filter_results(fp_val, fs_val, dp_val, ds_val, Fs_val, N_est, M_est, f_data, mag_data)
println("========================================")