clear;
N=input(' 滤波器阶数N =');
Wn=input(' 截止频率Wn = ');
Rp=input('通带波纹Rp = ');
[num,den]=cheby1(N,Rp,Wn,'s');
w=0:5*Wn;
h=freqs(num,den,w);
plot(w,20*log10(abs(h))),grid;
xlabel('Frequency, Hz'); ylabel('Gain, dB');
