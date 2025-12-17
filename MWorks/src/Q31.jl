using LinearAlgebra
using TyPlot  # 使用 Syslab 原生绘图库

# =========================================================================
#  第一部分：核心算法 (手动实现 residuez)
# =========================================================================

# 1. 求多项式根 (使用伴随矩阵法，数值稳定性高)
function get_roots_robust(coeffs)
    n = length(coeffs) - 1
    # 归一化，防止首项不为1
    if abs(coeffs[1]) < 1e-9; return ComplexF64[]; end
    c = coeffs ./ coeffs[1]
    
    # 构建伴随矩阵
    Cm = zeros(ComplexF64, n, n)
    for i in 1:n-1; Cm[i+1, i] = 1.0; end
    for i in 1:n; Cm[i, n] = -c[i+1]; end
    
    return eigvals(Cm)
end

# 2. 部分分式展开 (Residuez)
function my_residuez(b, a)
    # 计算直接项 k (长除法首项)
    k = 0.0
    if length(b) == length(a); k = b[1] / a[1]; end
    
    # 计算极点 p
    p = get_roots_robust(a)
    
    # 计算留数 r (使用留数定理导数法)
    r = zeros(ComplexF64, length(p))
    b_new = b .- k .* a # 去除直接项后的分子
    
    for i in 1:length(p)
        pi = p[i]
        # 分子值 Num(pi)
        num_val = sum(b_new[j] * (pi ^ (-(j-1))) for j in 1:length(b_new))
        # 分母导数值 Den'(pi) (等效于去除该极点后的连乘)
        den_val = 1.0 + 0.0im
        for j in 1:length(p)
            if i != j; den_val *= (1 - p[j]/pi); end
        end
        r[i] = num_val / den_val
    end
    return r, p, k
end

# =========================================================================
#  第二部分：主程序与结果输出
# =========================================================================

# 1. 输入题目给定的系数
# H(z) 分子: 2 + 5z^-1 + 1z^-2 - 3z^-3 + 4z^-4 + 6z^-5
b = [2.0, 5.0, 1.0, -3.0, 4.0, 6.0]
# H(z) 分母: 1 + 3z^-1 - 5z^-2 + 2z^-3 - 4z^-4 + 3z^-5
a = [1.0, 3.0, -5.0, 2.0, -4.0, 3.0]

println("正在进行 IIR 滤波器并联分解...\n")
r, p, k = my_residuez(b, a)

# 2. 打印计算结果 (作业答案)
println("=======================================================")
println("                分解结果 (可直接抄写)")
println("=======================================================")

# --- 输出并联结构 I 型 (复数形式) ---
println("\n[1] 并联结构 I 型 (Complex Form):")
println("    H(z) = k + Σ [ r_i / (1 - p_i z^-1) ]")
println("-------------------------------------------------------")
print("H(z) = ", round(real(k), digits=4))

for i in 1:length(p)
    # 格式化复数
    r_re = round(real(r[i]), digits=4); r_im = round(imag(r[i]), digits=4)
    r_str = abs(r_im)<1e-4 ? "$r_re" : "($r_re + $(r_im)j)"
    
    p_re = round(real(p[i]), digits=4); p_im = round(imag(p[i]), digits=4)
    p_str = abs(p_im)<1e-4 ? "$p_re" : "($p_re + $(p_im)j)"
    
    println("")
    print("       + [ $r_str ] / ( 1 - $p_str z^-1 )")
end
println("")

# --- 输出并联结构 II 型 (实数形式 - 二阶节合并) ---
println("\n\n[2] 并联结构 II 型 (Real Form):")
println("    说明: 将共轭复数极点合并为二阶节(SOS)，系数为实数")
println("-------------------------------------------------------")
print("H(z) = ", round(real(k), digits=4))

handled = falses(length(p)) # 标记已处理的极点

for i in 1:length(p)
    if handled[i]; continue; end
    
    # 判定实数极点
    if abs(imag(p[i])) < 1e-5
        val_r = real(r[i]); val_p = real(p[i])
        sign_p = val_p >= 0 ? "-" : "+" # 分母显示 1 - pz^-1
        
        println("")
        print("       + ", round(val_r, digits=4))
        print(" / ( 1 $sign_p ", round(abs(val_p), digits=4), " z^-1 )")
        handled[i] = true
    else
        # 寻找共轭对
        for j in (i+1):length(p)
            if !handled[j] && abs(real(p[i]) - real(p[j])) < 1e-5
                # 合并公式
                # 分子 b0 + b1*z^-1 = 2Re(r) - 2Re(r*p')z^-1
                b0 = 2 * real(r[i])
                b1 = -2 * real(r[i] * conj(p[i]))
                
                # 分母 1 + a1*z^-1 + a2*z^-2 = 1 - 2Re(p)z^-1 + |p|^2z^-2
                a1 = -2 * real(p[i])
                a2 = abs(p[i])^2
                
                # 格式化
                b0_s = round(b0, digits=4)
                b1_s = b1 >= 0 ? "+ $(round(b1,digits=4))" : "- $(round(abs(b1),digits=4))"
                a1_s = a1 >= 0 ? "+ $(round(a1,digits=4))" : "- $(round(abs(a1),digits=4))"
                a2_s = round(a2, digits=4)
                
                println("")
                print("       + ( $b0_s $b1_s z^-1 )") 
                print(" / ( 1 $a1_s z^-1 + $a2_s z^-2 )")
                
                handled[i] = true; handled[j] = true; break
            end
        end
    end
end
println("\n")

# =========================================================================
#  第三部分：绘图 (Pole-Zero Plot)
# =========================================================================

# 计算零点用于绘图
zeros_val = get_roots_robust(b)

figure("IIR Filter Analysis")

# 1. 绘制单位圆
theta = range(0, 2pi, length=200)
plot(cos.(theta), sin.(theta), "k--", linewidth=1, label="Unit Circle")
hold("on")

# 2. 绘制零点 (Zeros) - 蓝色圆圈
plot(real(zeros_val), imag(zeros_val), "bo", markersize=8, label="Zeros")

# 3. 绘制极点 (Poles) - 红色叉号
plot(real(p), imag(p), "rx", markersize=10, label="Poles")

# 4. 图表修饰
title("Pole-Zero Plot of IIR Filter")
xlabel("Real Axis")
ylabel("Imaginary Axis")
grid("on")
axis("equal") # 保证圆是圆的
legend("on")

hold("off")