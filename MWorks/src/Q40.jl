using TyPlot
using Random

# ===========================
# 1. 参数设置
# ===========================
n = 0:50                # 时间序列
M = 60                  # 平均次数
noise_amp = 0.6         # 噪声幅度

# ===========================
# 2. 生成原始信号 s[n]
# ===========================
s = 3 .* n .* (0.8) .^ n

# ===========================
# 3. 总体均值滤波模拟
# ===========================
x_accumulated = zeros(length(n)) 
x_single = zeros(length(n))      

println("正在进行计算...")

for i in 1:M
    # 生成噪声: 范围 [-0.6, 0.6]
    d = noise_amp .* (2 .* rand(length(n)) .- 1)
    
    # 受干扰信号
    x_current = s .+ d
    
    # 【修正点】记录第1次结果
    if i == 1
        # 使用 [:] 进行原地更新，确保修改的是外部变量
        x_single[:] = x_current
    end
    
    # 累加信号 (使用 global 也是一种方法，但对于数组推荐用原地更新)
    global x_accumulated += x_current
end

# 计算平均值
x_avg = x_accumulated ./ M

# ===========================
# 4. 绘图
# ===========================
figure()

# --- 子图1: 单次噪声序列 ---
subplot(2, 1, 1)
plot(n, x_single, "b.-") 
title("单次受干扰信号 (带噪声)")
xlabel("n")
ylabel("幅度")
grid("on")
xlim([0, 50])

# --- 子图2: 总体平均序列 ---
subplot(2, 1, 2)
plot(n, x_avg, "r.-", linewidth=1.5)
hold("on")
plot(n, s, "k--", linewidth=1.0) # 原始信号参考
title("60次检测结果的总体平均序列")
xlabel("n")
ylabel("幅度")
legend(["总体平均信号", "原始信号(参考)"])
grid("on")
xlim([0, 50])

println("绘图完成。现在你应该能看到上面的图有明显的噪声毛刺了。")