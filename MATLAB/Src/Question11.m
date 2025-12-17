clc;
A=[1 3.2 1.5 -0.8 1.4];
B=[1 -0.2 0.5];
[H,T]=impz(B,A,30);
disp(H);
stem(T,H);
