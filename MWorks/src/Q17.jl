using TyPlot
using Statistics # 导入统计包以使用 median 函数

# 1. 信号生成参数
N = 40                  # 序列长度
n = 0:N-1               # 时间索引 (0 到 39)

# 2. 生成原始信号 s[n]
# 公式: s[n] = 3 * n * (0.8)^n
# 注意: Julia 中 .^ 用于逐元素幂运算
s = 3 .* n .* (0.8).^n

# 3. 生成加性噪声 d[n]
# "幅度 0.6" 通常指均匀分布在 [-0.6, 0.6] 之间的随机噪声
# rand(N) 生成 [0, 1) 的随机数
# (rand(N) .- 0.5) .* 2 生成 [-1, 1)
# 乘以 0.6 得到 [-0.6, 0.6)
noise_amp = 0.6
d = (rand(N) .- 0.5) .* 2 .* noise_amp

# 生成受干扰序列 x[n]
x = s .+ d

# 4. 实现长度为 5 的中值滤波器
# 窗口长度 L = 5，意味着取 n-2, n-1, n, n+1, n+2 的中位数
L = 5
half_L = floor(Int, L / 2) # 半窗口长度 = 2
y = zeros(N) # 初始化输出序列

for i in 1:N
    # 确定当前窗口的起始和结束索引 (处理边界情况)
    # Julia 索引从 1 开始，对应 n = i-1
    # 窗口范围: [max(1, i - half_L), min(N, i + half_L)]
    start_idx = max(1, i - half_L)
    end_idx = min(N, i + half_L)
    
    # 提取窗口内的数据
    window_data = x[start_idx:end_idx]
    
    # 计算中值并赋值
    y[i] = median(window_data)
end

# 5. 绘图
figure("Median Filter Analysis")

# 子图 1: 受干扰的序列 x[n]
subplot(2, 1, 1)
stem(n, x, "b-o", label="受干扰序列")
# 可选: 绘制原始纯净信号作为参考 (虚线)
plot(n, s, "g--", label="原始信号") 
title("受干扰序列 (原始信号 + 噪声)")
ylabel("幅度")
legend()
grid(true)

# 子图 2: 中值滤波器输出 y[n]
subplot(2, 1, 2)
stem(n, y, "r-o", label="滤波后输出")
title("5点中值滤波器输出")
xlabel("样本 n")
ylabel("幅度")
legend()
grid(true)