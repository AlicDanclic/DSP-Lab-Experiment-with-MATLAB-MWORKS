% Design of IIR Butterworth Bandpass Filter
Wp =[0.4 0.6];
Ws = [0.3 0.7];
Rp = 0.6;
Rs = 35;
[N,Wn] = buttord(Wp, Ws, Rp, Rs);
[b,a] = butter(N,Wn);
[h,omega] = freqz(b,a,256);
plot (omega/pi,abs(h));grid;
xlabel('\omega/\pi'); ylabel('Gain');
title('IIR Butterworth Bandpass Filter');
disp(N);
disp(Wn);
