clear;
N=input('Type in the order N = ');
Wn=input('Type in the 3-dB cutoff frequency Wn = '); %模拟频率
[num,den]=butter(N,Wn,'s');
w=0:2*Wn;
h=freqs(num,den,w);
plot(w,20*log(abs(h))),grid;