% Program 10_5
% Lowpass Filter Design Using the Kaiser Window
clc
%delta_s=50;   %阻带衰减
delta_p=0.04;    %通带衰减
%alpha_s=10^(-delta_s/20);
alpha_p=1-10^(-delta_p/20);
fpts = [0.45 0.55];  %截止频率1，截止频率2
mag = [0 1];       %截止频率1对应的幅度 截止频率2对应的幅度
dev = [alpha_p alpha_p]; %通带衰减  阻带衰减（线性）
[N,Wn,beta,ftype]=kaiserord(fpts,mag,dev);
b = fir1(N,Wn,'high',kaiser(N+1,beta));
[h,omega] = freqz(b,1,512);
plot(omega/pi,20*log10(abs(h)));grid;
xlabel('\omega/\pi'); ylabel('Gain, dB');
