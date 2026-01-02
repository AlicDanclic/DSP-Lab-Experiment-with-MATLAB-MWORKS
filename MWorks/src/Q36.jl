using TyPlot
using LinearAlgebra

# DFT（沿用你文件里的实现）
function my_dft(x)
    N = length(x)
    X = zeros(ComplexF64, N)
    for k in 0:N-1
        s = 0.0 + 0.0im
        for n in 0:N-1
            s += x[n+1] * exp(-2im*pi*k*n/N)
        end
        X[k+1] = s
    end
    X
end

N = 64
n = collect(0:N-1)
x = sin.(25*pi*n/64)

Xk = my_dft(x)

# 幅度（不调用 abs）：|X| = sqrt(Re^2 + Im^2)
mag = sqrt.(real.(Xk).^2 .+ imag.(Xk).^2)

k = collect(0:N-1)

figure("Q36: 64-point DFT of x[n]")

subplot(3,1,1)
stem(k, real.(Xk), "b-", filled=true, markersize=4)
title("Re{X[k]}")
xlabel("k"); ylabel("Real")
grid("on")

subplot(3,1,2)
stem(k, imag.(Xk), "b-", filled=true, markersize=4)
title("Im{X[k]}")
xlabel("k"); ylabel("Imag")
grid("on")

subplot(3,1,3)
stem(k, mag, "b-", filled=true, markersize=4)
title("Magnitude |X[k]|  (computed by sqrt(Re^2+Im^2))")
xlabel("k"); ylabel("Magnitude")
grid("on")

# 标记真实“半个谱线”的位置：k0=12.5 和 N-k0=51.5
hold("on")
plot([12.5, 12.5], [0, maximum(mag)], "r--", linewidth=1, label="k=12.5")
plot([51.5, 51.5], [0, maximum(mag)], "r--", linewidth=1, label="k=51.5")
legend("on")
hold("off")