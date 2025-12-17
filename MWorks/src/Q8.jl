import Pkg

# ---------------------------------------------------------
# 依赖检查
# ---------------------------------------------------------
# 检查并安装 DSP 库
try
    using DSP
catch
    println("正在安装 DSP 包...")
    Pkg.add("DSP")
    using DSP
end

# 加载 TyPlot (MWorks 原生绘图库)
using TyPlot

# ---------------------------------------------------------
# 1. 定义滤波器性能指标
# ---------------------------------------------------------
# 归一化频率 (ω/π)
# 注意：DSP.jl 的 buttord 函数要求带通频率范围必须是 Tuple (元组)，不能是 Vector [向量]
wp = (0.4, 0.6)      # 通带截止频率 (使用圆括号)
ws = (0.3, 0.7)      # 阻带截止频率 (使用圆括号)
Rp = 0.6             # 通带最大波纹 (dB)
Rs = 35              # 阻带最小衰减 (dB)

println("正在计算滤波器参数...")

# ---------------------------------------------------------
# 2. 计算巴特沃兹滤波器参数
# ---------------------------------------------------------
# buttord(wp::Tuple, ws::Tuple, Rp, Rs)
N, Wn = buttord(wp, ws, Rp, Rs)

println("----------- 设计结果 -----------")
println("滤波器类型: Butterworth IIR Bandpass")
println("计算得到的阶数 N: ", N)
println("3dB 截止频率 Wn: ", Wn)

# ---------------------------------------------------------
# 3. 设计滤波器
# ---------------------------------------------------------
# Wn 返回的也是 Tuple，直接通过索引访问
responsetype = Bandpass(Wn[1], Wn[2])
designmethod = Butterworth(N)
filter_system = digitalfilter(responsetype, designmethod)

# ---------------------------------------------------------
# 4. 分析频率响应
# ---------------------------------------------------------
num_points = 1024
w_vals = range(0, π, length=num_points)

# 计算频率响应
H = freqz(filter_system, w_vals)
mag_db = 20 * log10.(abs.(H))

# x轴数据 (归一化频率)
x_axis = w_vals / π

# ---------------------------------------------------------
# 5. 使用 TyPlot 绘制增益响应曲线
# ---------------------------------------------------------
# 创建图形窗口
figure("Butterworth IIR Bandpass Filter")

# 绘制主响应曲线
plot(x_axis, mag_db, "b-", label="Filter Response", linewidth=2)

# 添加网格 (修改为 "on")
grid("on")

# 设置标题和轴标签
title("Butterworth IIR Bandpass Filter Gain Response")
xlabel("Normalized Frequency (×π rad/sample)")
ylabel("Magnitude (dB)")

# 保持当前图形以便添加辅助线 (修改为 "on")
hold("on")

# 添加指标辅助线 (通带 -0.6dB)
plot([0.4, 0.6], [-Rp, -Rp], "g--", label="Passband Spec (-0.6dB)", linewidth=1.5)

# 添加指标辅助线 (阻带 -35dB)
# 分两段画，避免中间连线
plot([0.0, 0.3], [-Rs, -Rs], "r--", label="Stopband Spec (-35dB)", linewidth=1.5)
plot([0.7, 1.0], [-Rs, -Rs], "r--", linewidth=1.5)

# 设置Y轴范围
ylim(-80, 5)

# 显示图例
legend()

println("绘图完成。")