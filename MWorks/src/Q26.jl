using TyPlot
using FFTW # 必须引入 FFTW 库来计算 DFT/IDFT

# 1. 定义已知序列
x1 = [2.2, 3.0, -1.5, 4.2, -1.8]
x2 = [0.8, -1.0, 1.6, 0.8]

# 2. 确定卷积后的长度 L
# 线性卷积的长度 = N1 + N2 - 1
N1 = length(x1)
N2 = length(x2)
L = N1 + N2 - 1

println("序列 x1 长度: $N1")
println("序列 x2 长度: $N2")
println("线性卷积所需最小 DFT 长度 L: $L")

# 3. 对序列进行补零 (Zero Padding)
# 将序列延长到长度 L，不足的部分补 0
x1_padded = [x1; zeros(L - N1)]
x2_padded = [x2; zeros(L - N2)]

# 4. 执行 DFT (FFT 算法)
X1_k = fft(x1_padded)
X2_k = fft(x2_padded)

# 5. 在频域相乘
# 卷积定理: 时域卷积 <-> 频域乘积
X_k = X1_k .* X2_k

# 6. 执行 IDFT (逆变换) 并取实部
# ifft 得到的结果通常包含极小的虚部误差，需要取 real()
x_conv = real(ifft(X_k))

# 7. 打印结果
println("基于 DFT 计算的卷积结果 x[n]:")
println(x_conv)

# 8. 绘图
figure("Convolution via DFT")

# 定义时间轴 n
n = 0:(L-1)

# 使用 stem 绘制离散序列
stem(n, x_conv)

title("基于 DFT 计算的线性卷积 x[n]")
xlabel("n (样本)")
ylabel("幅度")
grid("on")