clc;
clear;
x1=[2.2 3 -1.5 4.2 -1.8];
x2=[0.8 -1 1.6 0.8];
n=size(x1,2);
m=size(x2,2);
X1=fft(x1,n+m-1);
X2=fft(x2,n+m-1);
X=X1.*X2;
x=ifft(X)
stem(x,'.');
