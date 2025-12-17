using DSP
using TyPlot

"""
计算并绘制离散系统的脉冲响应
系统函数: H(z) = (1 - 0.2z^-1 + 0.5z^-2) / (1 + 3.2z^-1 + 1.5z^-2 - 0.8z^-3 + 1.4z^-4)
"""
function solve_impulse_response()
    # 1. 定义系统系数
    # 分子系数 b (对应 z^-0, z^-1, z^-2)
    b = [1.0, -0.2, 0.5]
    
    # 分母系数 a (对应 z^-0, z^-1, z^-2, z^-3, z^-4)
    a = [1.0, 3.2, 1.5, -0.8, 1.4]
    
    # 2. 生成输入信号：单位脉冲 (Unit Impulse)
    # 我们需要前 30 个样本
    N = 30
    delta = zeros(Float64, N)
    delta[1] = 1.0  # 在 n=0 处为 1 (Julia 索引从 1 开始，所以对应数组第 1 个元素)

    # 3. 计算脉冲响应
    # 使用 DSP.filt 函数对单位脉冲进行滤波，得到的输出即为脉冲响应 h[n]
    h = DSP.filt(b, a, delta)

    # 4. 打印数值结果
    println("=== 脉冲响应 h[n] 前 30 个样本值 ===")
    for i in 1:N
        # n 从 0 开始，数组索引 i 从 1 开始
        println("n = $(i-1):  $(h[i])")
    end

    # 5. 绘制图形 (离散序列通常使用杆图 stem，这里用带标记的线图模拟)
    TyPlot.clf()
    
    # 创建时间轴 n = 0 到 29
    n_axis = 0:(N-1)
    
    # 绘制
    TyPlot.plot(n_axis, h, "bo-", linewidth=1, markersize=5, label="h[n]")
    
    # 设置图表属性
    TyPlot.title("系统脉冲响应 h[n] (Impulse Response)")
    TyPlot.xlabel("样本序号 n")
    TyPlot.ylabel("幅值")
    TyPlot.grid(true)
    TyPlot.legend()
    
    println("\n计算与绘图完成。")
end

# 运行函数
solve_impulse_response()