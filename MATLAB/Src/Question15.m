num=[1,-0.2,0.5];
den=[1,3.2,1.5,-0.8,1.4];
x=[1 zeros(1,20-1)];
y=filter(num, den, x);
disp('Coefficients of the power series expansion');
disp(y)
stem(y)
