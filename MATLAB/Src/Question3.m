% Partial-Fraction Expansion of Rational z-Transform
num = [0 0 1 -0.2 0.5];
den = [1 3.2 1.5 -0.8 1.4];
[r,p,k] = residuez(num,den);
disp('Residues');disp(r')
disp('Poles');disp(p')
disp('Constants');disp(k)