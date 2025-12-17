n=0:15;x=cos(pi*n*0.5);
X=fft(x);
subplot(2,1,1);
stem(n,X,'.');
title('Magnitude of DFT')
xlabel('n'); ylabel('Magnitude')
%circulating DTFT
k=0:499;
w = pi/500*k;
X1=x*(exp(-1i*pi/500)).^(n'*k);
magX=abs(X1);
subplot(2,1,2);
plot(w/pi,magX);title('幅度响应');grid;
ylabel('幅度');xlabel('以\pi为单位的频率');
