N = 4;
Wn = 1;
[num,den] = butter(N,Wn,'s');
[h,w] = freqs(num,den);
plot (w,20*log10(abs(h)));
xlabel('Frequency, Hz'); ylabel('Gain, dB');
title('The 4th-order IIR Butterworth Lowpass Filter ')
grid on