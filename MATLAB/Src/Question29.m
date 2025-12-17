h = ones(1,6)/6; 
[H,w] = freqz(h, 1, 256);
subplot(211)
plot(w/pi, abs(H));
ylabel('Magnitude'); xlabel('\omega/\pi');
subplot(212)
ph = angle(H)*180/pi;
plot(w/pi,angle(H)*180/pi);
ylabel('Phase, degrees');xlabel('\omega/\pi');
