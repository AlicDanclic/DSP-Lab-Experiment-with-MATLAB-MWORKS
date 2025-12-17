xn=ones(10,1);
Xk=fft(xn,16);Xkf=abs(Xk);Xkp=angle(Xk);
subplot(211);
stem(0:15,Xkf,'filled');
xlabel('Time index/n');ylabel('Magnitude');
subplot(212);
stem(0:15,Xkp,'filled');
xlabel('Time index/n');ylabel('Phase')
