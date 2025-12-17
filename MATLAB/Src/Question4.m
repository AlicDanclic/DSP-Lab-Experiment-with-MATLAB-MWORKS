% a4
disp('prewapping is done,and T=2');
Wp = tan(0.75*pi/2);
Ws = tan(0.5*pi/2);
Rp = 0.5;
Rs = 43;
[N,Wn] = cheb1ord(Ws,Wp,Rp,Rs,'s');
[b,a] = cheby1(N,Rp,Wn,'s');
[bt,at]=lp2hp(b,a,Wp);
[num,den]=bilinear(bt,at,0.5);
[h,omega] = freqz(num,den);
plot (omega/pi,20*log10(abs(h)));grid;
xlabel('\omega/\pi'); ylabel('Gain');
title('Type I Chebyshev Highpass Filter');

clear;%预畸变
Rp=0.5;
Rs=43;
Wp=0.75;
Ws=0.35;
[N,Wp]=cheb1ord(Wp,Ws,Rp,Rs);
[num,den]=cheby1(N,Rp,Wp,'high');
w=0:pi/1024:pi;
h=freqz(num,den,w);
subplot(2,1,1);
plot(w/pi,abs(h)),grid;title('Amplitude in linear scale')
subplot(2,1,2);
plot(w/pi,20*log10(abs(h))),grid;
title('Amplitude in log scale')