%Q1_solution
%ellip(N,Ap,Ast,Wp)
%N--->The order of the filter
%Ap-->ripple in the passband
%Ast->a stopband Rs dB down from the peak value in the passband
%Wp-->the passband width
[be,ae]=ellip(3,0.6,32,0.4);
hellip=dfilt.df2(be,ae);
f=0:0.001:0.4;
g=grpdelay(hellip,f,2);
g1=max(g)-g;
[b,a,tau]=iirgrpdelay(6,f,[0 0.4],g1);
hallpass=dfilt.df2(b,a);
hoverall=cascade(hallpass,hellip);
hFVT=fvtool([hellip,hoverall]);
set(hFVT,'Filter',[hellip,hoverall]);
legend(hFVT,'Lowpass Elliptic filter','Compensated filter');

clear;
[num1,den1]=ellip(3,0.6,32,0.4);
[GdH,w]=grpdelay(num1,den1,512);
plot(w/pi,GdH); grid
xlabel('\omega/\pi'); ylabel('Group delay, samples');
F=0:0.001:0.4;
g=grpdelay(num1,den1,F,2);   % Equalize the passband
Gd=max(g)-g;
% Design the allpass delay equalizer
[num2,den2]=iirgrpdelay(6,F,[0,0.4],Gd);
[GdA,w] = grpdelay(num2,den2,512);
hold on;
plot(w/pi,GdH+GdA,'r'); 
legend('Original Filter','Compensated filter');
