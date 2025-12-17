clear
k=500;                      % number of frequency samples is 500
num=[1 -.2 .5 2 -.6];           %Numerator coefficients
den=[1 3.2 1.5 -.8 1.4];         %Denominator coefficients
w=0:pi/(k-1):pi;                 
h=freqz(num,den,w);             %Compute the frequency response
subplot(1,2,1)
plot(w/pi,abs(h))
title('Magnitude Spectrum')
xlabel('\omega/\pi');ylabel('Magnitude')
subplot(1,2,2)
plot(w/pi,unwrap(angle(h)))               %unwrapped phase function
title('Phase Spectrum')
xlabel('\omega/\pi');ylabel('Phase,radians')
[sos,g]=tf2sos(num,den)       
