% Program 10_2
% Design of Equiripple Linear-Phase FIR Filters
%
format long
fedge = [1500 1800];
FT = 5000;
mval = [1 0];
dev =[0.015 0.021];
[N,fpts,mag,wt] = remezord(fedge,mval,dev,FT);
b = remez(N,fpts,mag,wt);
disp('FIR Filter Coefficients'); disp(b)
[h,w] = freqz(b,1,256);
subplot(2,1,1);
plot(w/pi,20*log10(abs(h)));grid;
xlabel('\omega/\pi'); ylabel('Gain, dB');
subplot(2,1,2);
plot(w/pi,20*log10(abs(h)));grid;
axis([0 0.4 -0.7 0.7]);
xlabel('\omega/\pi'); ylabel('Gain, dB');
title('detail in passband')
