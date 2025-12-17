% Program 2_5
% Illustration of Median Filtering
%
N = 5;
R = 40; 
b = 6/5*(rand(1,R)-0.5); % Generate impulse noise
m = 0:R-1;
s = 3*m.*(0.8.^m); % Generate signal
x = s + b; % Impulse noise corrupted signal
y = medfilt1(x,N); % Median filtering
subplot(2,1,1)
stem(m,x);axis([0 50 -1 8]);grid on;
xlabel('n');ylabel('Amplitude');
title('Impulse Noise Corrupted Signal');
subplot(2,1,2)
stem(m,y);grid on;
xlabel('n');ylabel('Amplitude');
title('Output of Median Filter');
