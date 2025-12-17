N=64;
n=0:1:63;
x=sin(25*pi*n/N);
%k=512;
%w = 0:pi/(k-1):pi;
%h = freqz(x, 1, w);
%subplot(211);
%plot(w/pi,abs(h));grid
%title('Magnitude Spectrum')
%xlabel('\omega/\pi'); ylabel('Magnitude')

X=fft(x,64);
%subplot(212)
stem(n,X,'.');grid;