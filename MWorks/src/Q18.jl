using TyPlot
using FFTW # 需要 FFTW 包来进行快速傅里叶变换

# 1. 定义参数和频域序列 X[k]
N = 16
k = 0:N-1

# 根据题目定义: X[k] = k/16, for 0 <= k <= 15
X_k = k ./ 16.0

# 注意: 这里的 X[k] 是纯实数且不满足共轭对称性 (即 X[k] != conj(X[N-k]))
# 因此，反变换后的 x[n] 将是复数。

# 2. 计算 IDFT 得到 x[n]
# ifft 计算公式为: x[n] = (1/N) * sum(X[k] * exp(j*2*pi*n*k/N))
x_n = ifft(X_k)

# 3. 提取实部和虚部
x_real = real(x_n)
x_imag = imag(x_n)

# 打印部分数值结果供参考
println("x[n] 的前 5 个值 (复数形式):")
for i in 1:5
    println("n=$(i-1): $(round(x_n[i], digits=4))")
end

# 4. 绘图
figure("Sequence x[n] Real and Imaginary Parts")

# 子图 1: 实部
subplot(2, 1, 1)
stem(k, x_real, "b-o", label="实部")
title("x[n] 的实部")
ylabel("Re{x[n]}")
grid(true)
legend()

# 子图 2: 虚部
subplot(2, 1, 2)
stem(k, x_imag, "r-d", label="虚部") # 使用不同的颜色和标记
title("x[n] 的虚部")
xlabel("n")
ylabel("Im{x[n]}")
grid(true)
legend()