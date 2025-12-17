using DSP
using TyPlot

"""
椭圆滤波器设计与群延时分析
任务：
1. 设计 5 阶椭圆低通滤波器 (0.35pi, 0.8dB ripple, 35dB attenuation)
2. 框架性展示 10 阶全通滤波器级联 (注：自动系数优化需额外算法支持)
3. 绘制群延时对比图
"""
function analyze_group_delay()
    println("=== 滤波器群延时分析 ===")

    # ==============================
    # 1. 设计椭圆低通滤波器 (Lowpass)
    # ==============================
    # 规格: N=5, Fp=0.35pi (归一化频率 0.35), Rp=0.8dB, Rs=35dB
    N_lp = 5
    Wn = 0.35  # 归一化频率 (1.0 = Nyquist, 即 pi)
    Rp = 0.8
    Rs = 35.0

    println("1. 设计椭圆低通滤波器: Order=$N_lp, Wn=$Wn, Rp=$Rp dB, Rs=$Rs dB")
    
    # 设计滤波器
    # DSP.Lowpass(Wn) 中的 Wn 是归一化频率 (0~1)
    lp_responsetype = DSP.Lowpass(Wn)
    lp_designmethod = DSP.Elliptic(N_lp, Rp, Rs)
    lp_filter = DSP.digitalfilter(lp_responsetype, lp_designmethod)

    # ==============================
    # 2. 构建全通滤波器 (Allpass Equalizer)
    # ==============================
    # 目标：设计 10 阶全通滤波器均衡通带群延时
    
    N_ap = 10
    println("2. 初始化全通滤波器结构 (Order=$N_ap)")
    
    # --- [占位符系数] ---
    # 目前设置为 [1.0, 0, ...] 代表直通
    a_ap = zeros(Float64, N_ap + 1)
    a_ap[1] = 1.0 
    
    # 全通滤波器的分子系数是分母系数的倒序
    b_ap = reverse(a_ap)
    
    # 创建全通滤波器对象
    ap_filter = DSP.PolynomialRatio(b_ap, a_ap)

    # ==============================
    # 3. 计算群延时 (Group Delay)
    # ==============================
    println("3. 计算群延时特性...")
    
    # 频率点数量
    n_points = 512
    
    # [修复] 显式生成频率向量 (0 到 pi)，避免 grpdelay 返回标量导致的 BoundsError
    # range 返回的是一个迭代器，collect 转换为数组
    w_rad = collect(range(0, π, length=n_points))
    
    # 计算低通滤波器的群延时
    # 当传入频率向量时，grpdelay 仅返回延时向量
    gd_lp = DSP.grpdelay(lp_filter, w_rad)
    
    # 计算全通滤波器的群延时
    gd_ap = DSP.grpdelay(ap_filter, w_rad)
    
    # 级联系统的群延时 = 低通延时 + 全通延时
    gd_total = gd_lp + gd_ap

    # [优化] 将频率归一化 (0~1) 以便与 Wn 比较和绘图
    w_norm = w_rad / π
    
    # 提取通带部分的索引用于绘图 (只看通带内部情况)
    # 因为 Wn 是 0.35 (归一化)，我们对比 w_norm
    passband_indices = w_norm .<= (Wn * 1.2) #稍微多画一点以便观察截止点
    
    # ==============================
    # 4. 绘制图形
    # ==============================
    println("4. 绘制群延时曲线...")
    TyPlot.clf() # 清除旧图
    
    # 绘制低通滤波器的群延时
    TyPlot.plot(w_norm[passband_indices], gd_lp[passband_indices], "b-", linewidth=2, label="Lowpass Filter Delay")
    
    # 绘制级联系统的群延时
    TyPlot.plot(w_norm[passband_indices], gd_total[passband_indices], "r--", linewidth=1.5, label="Cascaded (LP + AP) Delay")
    
    TyPlot.title("Group Delay in Passband")
    TyPlot.xlabel("Normalized Frequency (x pi rad/sample)")
    TyPlot.ylabel("Group Delay (samples)")
    TyPlot.grid(true)
    TyPlot.legend()
    
    # 标记截止频率
    # TyPlot.axvline(x=Wn, color="k", linestyle="--", label="Cutoff")

    println("完成！请查看绘图窗口。")
end

# 运行分析
analyze_group_delay()