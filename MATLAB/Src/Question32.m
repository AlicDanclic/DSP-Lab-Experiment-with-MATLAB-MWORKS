fpts = [0 0.32 0.35 0.65 0.68 1];
mval = [0.6 0.6 0.2 0.2 0.8 0.8];
b = fir2(100,fpts,mval);
[h,omega] = freqz(b,1,512);
plot(omega/pi,abs(h));grid;
xlabel('\omega/\pi'); ylabel('Magnitude');