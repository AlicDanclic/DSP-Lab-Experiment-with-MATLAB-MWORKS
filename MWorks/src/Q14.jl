using TyPlot
using FFTW

"""
DFT 分析演示
任务：
1. 生成长度为 10 的矩形序列
2. 计算 16 点 DFT (补零)
3. 绘制幅度和相位谱
"""
function dft_rectangular_window()
    println("=== 离散傅立叶变换 (DFT) 分析 ===")

    # ==============================
    # 1. 信号生成
    # ==============================
    # 长度为 10 的矩形序列 (Rectangular Sequence)
    # x[n] = 1, for 0 <= n < 10
    N_seq = 10
    x = ones(Float64, N_seq)
    
    println("原始序列长度: $N_seq")

    # ==============================
    # 2. 16 点 DFT 计算
    # ==============================
    N_fft = 16
    println("DFT 变换点数: $N_fft")

    # 构造补零后的序列
    # 将 x 放入长度为 16 的零数组的前 10 个位置
    x_padded = zeros(Float64, N_fft)
    x_padded[1:N_seq] = x

    # 计算 FFT
    # 注意：Julia 的 fft 函数位于 FFTW 包中
    X_k = fft(x_padded)

    # ==============================
    # 3. 计算幅度与相位
    # ==============================
    # 幅度 Spectrum Magnitude
    magnitude = abs.(X_k)

    # 相位 Spectrum Phase (弧度)
    phase = angle.(X_k)

    # ==============================
    # 4. 绘图
    # ==============================
    println("正在绘制幅度和相位谱...")
    TyPlot.clf() # 清除旧图

    # 频率轴 k = 0, 1, ..., 15
    k = 0:(N_fft-1)

    # --- 子图 1: 幅度谱 ---
    TyPlot.subplot(2, 1, 1)
    # 使用 "bo-" (蓝色圆点实线) 模拟离散频谱的杆状图效果
    TyPlot.plot(k, magnitude, "bo-", linewidth=1.5, markersize=5, label="|X[k]|")
    TyPlot.title("16点 DFT 幅度谱 (Magnitude)")
    TyPlot.xlabel("频率索引 k")
    TyPlot.ylabel("幅度")
    TyPlot.grid(true)
    # 调整 Y 轴范围以更美观地显示主瓣
    TyPlot.ylim(0, maximum(magnitude) * 1.1)

    # --- 子图 2: 相位谱 ---
    TyPlot.subplot(2, 1, 2)
    # 使用 "ro-" (红色圆点实线)
    TyPlot.plot(k, phase, "ro-", linewidth=1.5, markersize=5, label="∠X[k]")
    TyPlot.title("16点 DFT 相位谱 (Phase)")
    TyPlot.xlabel("频率索引 k")
    TyPlot.ylabel("相位 (rad)")
    TyPlot.grid(true)
    # 设置 Y 轴范围为 -pi 到 pi，方便观察
    TyPlot.ylim(-pi, pi)

    println("处理完成。")
end

# 运行分析函数
# 注意：如果提示 UndefVarError: `fft` not defined，请先运行 `import Pkg; Pkg.add("FFTW")`
try
    dft_rectangular_window()
catch e
    println("运行出错: ", e)
    println("请确保已安装 FFTW 包: import Pkg; Pkg.add(\"FFTW\")")
end