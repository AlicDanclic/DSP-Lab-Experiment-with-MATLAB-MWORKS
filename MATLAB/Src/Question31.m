num = [2 5 1 -3 4 6];
den = [1 3 -5 2 -4 3];
[r1,p1,k1] = residuez(num, den); %系统并联I 型分解
[r2,p2,k2] = residue(num, den); %系统并联II 型分解
disp('Parallel Form I')
disp('Residues are'); disp(r1);
disp('Poles are at'); disp(p1);
disp('Constant value'); disp(k1);
disp('Parallel Form II')
disp('Residues are'); disp(r2);
disp('Poles are at'); disp(p2);
disp('Constant value'); disp(k2);

%结果出来后，寻找其中的共轭复数对。制成rri=[ ] ppi=[ ] 按照题目中系数，构造如下
rri=[r1(3), r1(4)];
ppi=[p1(3), p1(4)];
[b1,a1]=residuez(rri, ppi, 0);
disp('Parallel Form I 复数极点顶展开式分子系数'); disp(b1);
disp('Parallel Form I 复数极点顶展开式分母系数'); disp(a1);

rr2=[r2(3), r2(4)];
pp2=[p2(3), p2(4)];
[b2,a2]=residue(rr2, pp2, 0);
disp('Parallel Form II 复数极点顶展开式分子系数'); disp(b2);
disp('Parallel Form II 复数极点顶展开式分母系数'); disp(a2);