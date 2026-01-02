using DSP
using TyPlot
using Optim
using Statistics
using LinearAlgebra

"""
第12题修正：设计 10 阶全通滤波器均衡群延时 (含优化算法)
"""
function design_and_analyze_group_delay()
    println("=== 第12题：群延时均衡设计 (优化中...) ===")

    # ==============================
    # 1. 设计椭圆低通滤波器
    # ==============================
    N_lp = 5
    Wn = 0.35  # 归一化频率
    Rp = 0.8
    Rs = 35.0
    
    lp_responsetype = DSP.Lowpass(Wn)
    lp_designmethod = DSP.Elliptic(N_lp, Rp, Rs)
    lp_filter = DSP.digitalfilter(lp_responsetype, lp_designmethod)

    # ==============================
    # 2. 优化全通均衡器参数
    # ==============================
    # 全通滤波器阶数 N_ap = 10 (5个二阶节)
    n_sections = 5 
    
    # 优化目标频率范围：只关注通带 (0 到 0.35pi)
    w_pass = range(0, stop=Wn*pi, length=100)
    tau_lp = grpdelay(lp_filter, w_pass)
    
    # 辅助函数：根据参数构建全通滤波器
    # 参数格式：[r1, theta1, r2, theta2, ...]
    function make_allpass(params)
        total_filter = nothing
        for i in 1:n_sections
            r = params[2*i - 1]
            theta = params[2*i]
            p = r * exp(im * theta)
            # 二阶全通节
            poles = [p, conj(p)]
            zeros = [1/conj(p), 1/p]
            # 为了保证实系数，增益通常处理为1，这里简化处理
            section = ZeroPoleGain(zeros, poles, 1.0)
            
            if total_filter === nothing
                total_filter = section
            else
                total_filter = total_filter * section
            end
        end
        return total_filter
    end

    # 代价函数：群延时标准差
    function cost_function(params)
        try
            ap = make_allpass(params)
            tau_ap = grpdelay(ap, w_pass)
            return std(tau_lp + tau_ap)
        catch
            return Inf
        end
    end

    # 初始猜测 & 边界
    initial_params = repeat([0.8, 0.2*pi], n_sections)
    # 稍微扰动一下初始值避免对称性陷阱
    for i in 1:length(initial_params); initial_params[i] += 0.01*i; end
    
    lower = repeat([0.0, 0.0], n_sections)
    upper = repeat([0.99, pi], n_sections)

    # 执行优化
    res = optimize(cost_function, lower, upper, initial_params, Fminbox(BFGS()), Optim.Options(time_limit=15.0))
    best_params = Optim.minimizer(res)
    
    # ==============================
    # 3. 结果绘图
    # ==============================
    ap_filter = make_allpass(best_params)
    
    # 绘图频率轴
    w_plot = range(0, stop=Wn*pi, length=300)
    tau_lp_plot = grpdelay(lp_filter, w_plot)
    tau_ap_plot = grpdelay(ap_filter, w_plot)
    tau_total = tau_lp_plot + tau_ap_plot
    
    TyPlot.clf()
    w_norm = w_plot ./ pi
    
    TyPlot.plot(w_norm, tau_lp_plot, "b--", linewidth=1, label="Original LP Delay")
    TyPlot.plot(w_norm, tau_total, "r-", linewidth=2, label="Equalized Total Delay")
    
    TyPlot.title("Group Delay Equalization (Order 10 Allpass)")
    TyPlot.xlabel("Normalized Frequency")
    TyPlot.ylabel("Group Delay (samples)")
    TyPlot.legend()
    TyPlot.grid(true)
    
    println("优化完成。通带群延时标准差从 $(round(std(tau_lp_plot), digits=2)) 降低到 $(round(std(tau_total), digits=2))")
end

design_and_analyze_group_delay()