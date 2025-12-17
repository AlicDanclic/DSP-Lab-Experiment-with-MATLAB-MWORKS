# 导入必要的库
using DSP
using TyPlot

"""
MWorks 椭圆低通滤波器设计工具
功能：
1. 通过键盘接收用户输入的滤波器参数
2. 计算椭圆滤波器的频率响应
3. 使用 TyPlot 绘制增益响应 (Bode Plot Magnitude)
"""
function design_and_plot_elliptic_filter()
    # ==============================
    # 1. 获取用户输入
    # ==============================
    println("=== 椭圆低通滤波器设计 (MWorks) ===")
    
    try
        # 强制刷新输出缓冲区，确保提示词在输入框之前显示
        flush(stdout)

        print("请输入滤波器阶数 (整数, 例如 4): ")
        input_str = readline()
        if isempty(input_str); println("输入为空"); return; end
        N = parse(Int, input_str)

        print("请输入采样频率 Fs (Hz, 例如 1000): ")
        input_str = readline()
        if isempty(input_str); println("输入为空"); return; end
        Fs = parse(Float64, input_str)

        print("请输入截止频率 Fc (Hz, 必须小于 Fs/2, 例如 200): ")
        input_str = readline()
        if isempty(input_str); println("输入为空"); return; end
        Fc = parse(Float64, input_str)

        print("请输入通带波纹 Rp (dB, 例如 1.0): ")
        input_str = readline()
        if isempty(input_str); println("输入为空"); return; end
        Rp = parse(Float64, input_str)

        print("请输入阻带衰减 Rs (dB, 例如 40.0): ")
        input_str = readline()
        if isempty(input_str); println("输入为空"); return; end
        Rs = parse(Float64, input_str)

        # 验证 Nyquist 频率限制
        if Fc >= Fs / 2
            println("错误：截止频率必须小于采样频率的一半 (Nyquist 频率)。")
            return
        end

        # ==============================
        # 2. 滤波器设计
        # ==============================
        println("\n正在计算滤波器系数...")
        
        # [修复] 使用 DSP.Lowpass 明确指定包名
        response_type = DSP.Lowpass(Fc; fs=Fs)
        
        # [修复] 使用 DSP.Elliptic 明确指定包名，解决 UndefVarError 问题
        design_method = DSP.Elliptic(N, Rp, Rs)
        
        # [修复] 使用 DSP.digitalfilter
        filter_object = DSP.digitalfilter(response_type, design_method)

        # ==============================
        # 3. 计算频率响应
        # ==============================
        # 生成频率点：从 0 到 Nyquist 频率，取 1024 个点
        freqs_range = range(0, Fs/2, length=1024)
        
        # [修复] 使用 DSP.freqz
        h = DSP.freqz(filter_object, freqs_range, Fs)
        
        # 计算幅值增益 (dB)
        magnitude_db = 20 * log10.(abs.(h))

        # ==============================
        # 4. 使用 TyPlot 绘图
        # ==============================
        println("正在绘制增益响应曲线...")

        # 清除当前图形 (如果有)
        TyPlot.clf()

        # 绘制曲线
        TyPlot.plot(freqs_range, magnitude_db, "b-", linewidth=1.5, label="Gain Response")
        
        # 设置标题和标签
        TyPlot.title("椭圆低通滤波器增益响应 (Elliptic Lowpass Filter)")
        TyPlot.xlabel("频率 (Hz)")
        TyPlot.ylabel("增益 (dB)")
        
        # 添加网格
        TyPlot.grid(true)
        
        # 添加图例
        TyPlot.legend()
        
        println("绘图完成！")
        
    catch e
        # 打印详细错误栈，方便调试
        println("运行出错: ", e)
        # showerror(stdout, e, catch_backtrace()) # 如果需要更详细的堆栈信息可以取消注释
        println("\n提示：请检查 DSP 包是否已安装 (import Pkg; Pkg.add(\"DSP\"))")
    end
end

# 运行主函数
design_and_plot_elliptic_filter()