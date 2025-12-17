clc;
clear;
x1=[2.2 3 -1.5 4.2 -1.8];   
x2=[0.8 -1 1.6 0.8];
x=conv(x1,x2)      %结果在主界面输出
stem(x,'filled');
grid on;
xlabel('Time index/n');ylabel('Amplitude');
title('The output convolution');
