format long
zr =[2.2 -1+i -1-i 1.4];
pr =[3.7+2*i 3.7-2*i -2.1-i -2.1+i ];
% Transpose zero and pole row vectors
z = zr'; p = pr';
k = 1;
[num, den] = zp2tf(z, p, k);
disp('Numerator polynomial coefficients'); disp(num);
disp('Denominator polynomial coefficients'); disp(den);