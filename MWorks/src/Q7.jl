# 引入包管理器
import Pkg

# 尝试引入 SymPy，如果报错则自动安装
try
    using SymPy
catch e
    println("正在安装 SymPy 包，请稍候...")
    Pkg.add("SymPy")
    using SymPy
end

# 定义符号 z
@vars z

# 定义 z^-1 (为了方便书写公式)
inv_z = z^(-1)

# 输入题目给定的系统函数 H(z)
# H(z) = 0.2 + 1/(1+3.2z^-1) + 0.6/(1-2.4z^-1) + 1.8/((1-2.4z^-1)^2)
H = 0.2 + 1/(1 + 3.2*inv_z) + 0.6/(1 - 2.4*inv_z) + 1.8/((1 - 2.4*inv_z)^2)

# 计算有理形式 (通分合并)
# together 函数会将多项式合并到同一个分母上
H_rational = together(H)

# 进一步简化表达式（展开分子分母）
# 这一步通常会将结果整理为 z 的正幂次形式
H_simplified = simplify(H_rational)

# 获取分子和分母
num = numer(H_simplified)
den = denom(H_simplified)

# 打印结果
println("----------- 计算结果 -----------")
println("原始表达式 H(z):")
println(H)
println("\n合并后的有理形式 H_rational(z):")
println(H_simplified)

println("\n----------- 分子与分母 -----------")
println("分子 N(z): ", 5.76*num)
println("分母 D(z): ", 5.76*den)

# 如果需要转换回 z^-1 的形式 (DSP 中常用)，可以上下同除以 z 的最高次幂
# 这里主要展示数学上的有理多项式形式