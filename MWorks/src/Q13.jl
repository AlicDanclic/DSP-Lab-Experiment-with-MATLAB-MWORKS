using TyPlot
using Random
using DSP

"""
4 点滑动平均滤波器演示
1. 生成原始信号 s[n] = 3 * n * (0.8)^n
2. 生成噪声信号 d[n] (幅度 0.6)
3. 合成受干扰信号 x[n] = s[n] + d[n]
4. 使用 4 点滑动平均滤波器进行处理
5. 绘制所有信号波形
"""
function moving_average_demo()
    # ==============================
    # 1. 信号生成
    # ==============================
    N = 40              # 序列长度
    n = 0:(N-1)         # 时间索引 n = 0, 1, ..., 39

    # 原始信号 s[n]
    s = 3 .* n .* (0.8).^n

    # 噪声信号 d[n]
    # 设定随机种子以保证结果可重复
    Random.seed!(42)
    # 生成幅度为 0.6 的随机噪声，假设为均匀分布 [-0.6, 0.6]
    d = 1.2 .* rand(N) .- 0.6

    # 受干扰信号 x[n]
    x = s + d

    # ==============================
    # 2. 滤波器设计与应用
    # ==============================
    # 4 点滑动平均滤波器
    # 差分方程: y[n] = 1/4 * (x[n] + x[n-1] + x[n-2] + x[n-3])
    # 系统函数: H(z) = 0.25 + 0.25z^-1 + 0.25z^-2 + 0.25z^-3
    
    M = 4
    b = ones(M) / M   # 分子系数 [0.25, 0.25, 0.25, 0.25]
    a = [1.0]         # 分母系数 [1.0]

    # 使用 DSP.filt 进行滤波
    y = DSP.filt(b, a, x)

    # ==============================
    # 3. 绘图
    # ==============================
    println("正在绘制波形...")
    TyPlot.clf() # 清除旧图

    # 使用 2x2 子图布局
    
    # 子图 1: 原始信号 s[n]
    TyPlot.subplot(2, 2, 1)
    TyPlot.plot(n, s, "b.-", label="s[n]") # 蓝色点线
    TyPlot.title("原始信号 s[n]")
    TyPlot.grid(true)
    # TyPlot.ylabel("幅值")

    # 子图 2: 噪声信号 d[n]
    TyPlot.subplot(2, 2, 2)
    TyPlot.plot(n, d, "g.-", label="d[n]") # 绿色点线
    TyPlot.title("噪声 d[n]")
    TyPlot.grid(true)

    # 子图 3: 受干扰信号 x[n]
    TyPlot.subplot(2, 2, 3)
    TyPlot.plot(n, x, "r.-", label="x[n]") # 红色点线
    TyPlot.title("受干扰信号 x[n]")
    TyPlot.xlabel("样本 n")
    TyPlot.ylabel("幅值")
    TyPlot.grid(true)

    # 子图 4: 滤波输出 y[n]
    TyPlot.subplot(2, 2, 4)
    TyPlot.plot(n, y, "k.-", label="y[n]", linewidth=1.5) # 黑色加粗点线
    TyPlot.title("4点滑动平均输出 y[n]")
    TyPlot.xlabel("样本 n")
    TyPlot.grid(true)

    # 调整整体布局标题 (如果支持)
    # TyPlot.suptitle("滑动平均滤波器效果对比")

    println("处理完成。")
    println("原始信号 s[n] (前5个点): ", s[1:5])
    println("滤波输出 y[n] (前5个点): ", y[1:5])
end

# 运行主函数
moving_average_demo()