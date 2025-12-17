% Program 2_4
% Signal Smoothing by a Moving-Average Filter
R = 40;
d = 6/5*(rand(1,R)-0.5);
m = 0:1:R-1;
s =3.*m.*0.8.^m;
x  = s + d;
subplot(211);
plot(m,d,'r-',m,s,'b:',m,x,'m--')
title('The sequence with noise');
ylabel('Amplitude')
legend('d[n]','s[n]','x[n]');
b = ones(4,1)/4;
y = fftfilt(b,x);
subplot(212);
plot(m,s,'r-',m,y,'b--') 
title('The original sequence & the output sequence');
legend('s[n]','y[n]');
ylabel('Amplitude')
