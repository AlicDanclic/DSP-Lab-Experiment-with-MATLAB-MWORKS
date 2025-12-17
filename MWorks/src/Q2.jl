import Pkg
# 检查并添加必要的包 (如果尚未安装)
required_packages = ["TyPlot"]
for pkg in required_packages
    if Base.find_package(pkg) === nothing
        println("正在安装 $pkg ...")
        Pkg.add(pkg)
    end
end

using TyPlot
using LinearAlgebra

println("=== 巴特沃兹模拟低通滤波器设计 ===")

# ==========================================
# 1. 键盘输入参数 (优化版)
# ==========================================
function get_valid_input(prompt::String, parse_func::Function)
    while true
        print(prompt)
        flush(stdout) # 确保提示文字立即显示
        
        try
            # 读取一行并去除首尾空格
            input_str = strip(readline())
            
            # 如果是空行（比如误触回车），跳过本次循环继续等待
            if isempty(input_str)
                continue
            end
            
            # 尝试解析
            value = parse_func(input_str)
            return value
        catch
            # 解析失败时提示，而不是直接使用默认值退出
            println("输入格式不正确，请重新输入。")
        end
    end
end

try
    # 获取 N
    global N = get_valid_input("请输入滤波器阶数 N (整数, 例如 4): ", x -> parse(Int, x))

    # 获取 fc
    global fc = get_valid_input("请输入 3dB 截止频率 fc (Hz, 例如 1000): ", x -> parse(Float64, x))

    println("\n正在设计: 阶数 N = $N, 截止频率 fc = $fc Hz")

catch e
    # 只有在发生严重系统错误（如中断）时才捕获
    println("\n发生未知错误，程序终止。")
    rethrow(e)
end

# ==========================================
# 2. 滤波器设计 (计算极点)
# ==========================================
# 巴特沃兹滤波器的极点分布在 S 平面的圆上
# 归一化角频率 omega_c = 2 * pi * fc
wc = 2 * pi * fc

# 极点公式: sk = wc * exp(j * pi * (2k + N - 1) / 2N)
# k = 1 到 N
poles = ComplexF64[]
for k in 1:N
    theta = pi * (2 * k + N - 1) / (2 * N)
    s_k = wc * exp(im * theta)
    push!(poles, s_k)
end

println("\n计算得到的极点 (S平面):")
for (i, p) in enumerate(poles)
    println("p$i: $(round(real(p), digits=2)) + $(round(imag(p), digits=2))j")
end

# ==========================================
# 3. 计算频率响应
# ==========================================
# 定义频率范围用于绘图：从 0 到 3倍截止频率
f_plot = range(0, stop=3*fc, length=1000)
w_plot = 2 * pi .* f_plot

# 巴特沃兹模拟低通滤波器的传输函数 H(s) = wc^N / product(s - pk)
# 为了计算增益 |H(jw)|，我们将 s 替换为 jw

# 分子 (对于低通滤波器，DC增益为1，分子等于极点的乘积的模，即 wc^N)
numerator = wc^N

gains = Float64[]

for w in w_plot
    s = im * w
    # 分母 = (s - p1)*(s - p2)*...
    denominator = prod(s .- poles)
    
    # H(s) = Num / Den
    H = numerator / denominator
    
    # 存入幅值
    push!(gains, abs(H))
end

# 转换为 dB
gains_db = 20 .* log10.(gains)

# ==========================================
# 4. 绘制增益响应
# ==========================================
println("\n正在绘制增益响应...")

# 绘制幅频响应
plot(f_plot, gains_db, "b-", linewidth=2, label="Gain Response")
hold("on")

# 标记 -3dB 截止频率点
plot([fc, fc], [-60, 5], "r--", label="Cutoff Frequency (-3dB)")
plot(fc, -3.01, "ro") # 在曲线上标记点 (理论值约为 -3.01 dB)

title("Butterworth Analog Lowpass Filter Response (N=$N, fc=$(Int(fc))Hz)")
xlabel("Frequency (Hz)")
ylabel("Magnitude (dB)")
grid("on")
legend()

# 设置Y轴范围使图表更清晰
ylim(-40, 5)

println("绘图完成。")