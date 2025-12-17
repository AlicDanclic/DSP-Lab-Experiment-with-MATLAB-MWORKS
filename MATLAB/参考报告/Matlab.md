# MATLAB 上机实验

1. 已知3阶椭圆IIR数字低通滤波器的性能指标为：通带截止频率0.4π，通带波纹为0.6dB，最小阻带衰减为32dB。设计一个6阶全通滤波器对其通带的群延时进行均衡。绘制低通滤波器和级联滤波器的群延时。

```matlab
% Q1_solution
% ellip(N,Ap,Ast,Wp)
% N-->The order of the filter
% Ap-->ripple in the passband
% Ast->a stopband Rs dB down from the peak value in the passband
% Wp-->the passband width

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
g=grpdelay(num1,den1,F,2); % Equalize the passband
Gd=max(g)-g;

% Design the allpass delay equalizer
[num2,den2]=iirgrpdelay(6,F,[0,0.4],Gd);
[GdA,w] = grpdelay(num2,den2,512);
hold on;
plot(w/pi,GdH+GdA,'r');
legend('Original Filter','Compensated filter');
```

2. 设计巴特沃兹模拟低通滤波器，其滤波器的阶数和3-dB截止频率由键盘输入，程序能根据输入的参数，绘制滤波器的增益响应。

```matlab
clear;
N=input('Type in the order N = ');
Wn=input('Type in the 3-dB cutoff frequency Wn = '); %模拟频率
[num,den]=butter(N,Wn,'s');
w=0:2*Wn;
h=freqs(num,den,w);
plot(w,20*log(abs(h))),grid;
```

3. 已知系统的系统函数为：

$$H(z) = \frac{z^2 - 0.2z + 0.5}{z^4 + 3.2z^3 + 1.5z^2 - 0.8z + 1.4}$$

用MATLAB进行部分分式展开，并写出展开后的表达式。

```matlab
% Partial-Fraction Expansion of Rational z-Transform
num = [0 0 1 -0.2 0.5];
den = [1 3.2 1.5 -0.8 1.4];
[r,p,k] = residuez(num,den);
disp('Residues');disp(r')
disp('Poles');disp(p')
disp('Constants');disp(k)
```

4. 设计切比雪夫I型IIR数字高通滤波器，其性能指标为：通带波纹0.5dB，最小阻带衰减43dB，通带和阻带边缘频率0.75π和0.5π，绘制所设计的滤波器增益响应。

```matlab
% a4
disp('prewapping is done,and T=2');
Wp = tan(0.75*pi/2);
Ws = tan(0.5*pi/2);
Rp = 0.5;
Rs = 43;
[N,Wn] = cheb1ord(Ws,Wp,Rp,Rs,'s');
[b,a] = cheby1(N,Rp,Wn,'s');
[bt,at]=lp2hp(b,a,Wp);
[num,den]=bilinear(bt,at,0.5);
[h,omega] = freqz(num,den);
plot (omega/pi,20*log10(abs(h)));grid;
xlabel('\omega/\pi'); ylabel('Gain');
title('Type I Chebyshev Highpass Filter');

clear;%预畸变
Rp=0.5;
Rs=43;
Wp=0.75;
Ws=0.35;
[N,Wp]=cheb1ord(Wp,Ws,Rp,Rs);
[num,den]=cheby1(N,Rp,Wp,'high');
w=0:pi/1024:pi;
h=freqz(num,den,w);
subplot(2,1,1);
plot(w/pi,abs(h)),grid;title('Amplitude in linear scale')
subplot(2,1,2);
plot(w/pi,20*log10(abs(h))),grid;
title('Amplitude in log scale')
```

4. 已知复指数序列为：$x[n] = 0.2e^{(0.4+j0.5)n}$，绘制30点该序列的实部和虚部。

```matlab
n=0:29;
x=0.2*exp((0.4+1i*0.5)*n);
subplot(211);
stem(n,real(x));
xlabel('n');ylabel('real part');
grid on;
subplot(212);
stem(n,imag(x));
xlabel('n');ylabel('imag part');
grid on;
```

5. 设计切比雪夫I型模拟低通滤波器，其滤波器的阶数，3-dB截止频率和通带的波纹由键盘输入，程序能根据输入的参数，绘制滤波器的增益响应。

```matlab
clear;
N=input(' 滤波器阶数N =');
Wn=input(' 截止频率Wn = ');
Rp=input('通带波纹Rp = ');
[num,den]=cheby1(N,Rp,Wn,'s');
w=0:5*Wn;
h=freqs(num,den,w);
plot(w,20*log10(abs(h))),grid;
xlabel('Frequency, Hz'); ylabel('Gain, dB');
```

7. 已知系统的系统函数为：

$$H(z) = 0.2 + \frac{1}{1+3.2z^{-1}} + \frac{0.6}{1-2.4z^{-1}} + \frac{1.8}{(1-2.4z^{-1})^2}$$

用MATLAB求系统z变换的有理形式，并写出有理形式的表达式。

```matlab
r=[1 0.6 1.8];
p=[-3.2 2.4 2.4];
k=0.2;
[num, den] = residuez(r,p,k)
```

8. 设计巴特沃兹IIR数字带通滤波器，其性能指标为：归一化通带截止频率为[0.4, 0.6]，归一化阻带截止频率为[0.3, 0.7]，通带波纹为0.6dB，最小阻带衰减为35dB。绘制所设计的滤波器增益响应。

```matlab
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
```

8. 已知指数序列为：$x[n] = 2(0.9)^n$，绘制24点该序列。

```matlab
n=0:23;
x=2*0.9.^n;
stem(n,x,'.');
grid on;
ylabel('Amplitude');
xlabel('Time index');
```

9. 设计椭圆模拟低通滤波器，其滤波器的阶数，3-dB截止频率，通带的波纹和阻带衰减由键盘输入，程序能根据输入的参数，绘制滤波器的增益响应。

```matlab
clear;
N=input('Type in the order N = ');
Wn=input('Type in the 3-dB cutoff frequency Wn = ');
Rp=input('Type in the the passband ripple Rp = ');
Rs=input('Type in the the minimum stopband attenuation Rs = ');
[num,den]=ellip(N,Rp,Rs,Wn,'s');
w=0:5*Wn;
h=freqs(num,den,w);
plot(w,20*log10(abs(h))),grid;
xlabel('Frequency, Hz'); ylabel('Gain, dB');
```

11. 已知系统的系统函数为：

$$H(z) = \frac{z^2 - 0.2z + 0.5}{z^4 + 3.2z^3 + 1.5z^2 - 0.8z + 1.4}$$

用MATLAB的impz函数求h[n]的前30个样本值。

```matlab
clc;
A=[1 3.2 1.5 -0.8 1.4];
B=[1 -0.2 0.5];
[H,T]=impz(B,A,30);
disp(H);
stem(T,H);
```

12. 已知5阶椭圆IIR数字低通滤波器的性能指标为：通带截止频率0.35π，通带波纹为0.8dB，最小阻带衰减为35dB。设计一个10阶全通滤波器对其通带的群延时进行均衡。绘制低通滤波器和级联滤波器的群延时。

```matlab
%ellip(N,Ap,Ast,Wp)
%N-->The order of the filter
%Ap-->ripple in the passband
%Ast->a stopband Rs dB down from the peak value in the passband
%Wp-->the passband width

[be,ae]=ellip(5,0.8,35,0.35);
hellip=dfilt.df2(be,ae);
f=0:0.001:0.4;
g=grpdelay(hellip,f,2);
g1=max(g)-g;
[b,a,tau]=iirgrpdelay(10,f,[0 0.4],g1);
%the first parameter above is the order of the allpass
hallpass=dfilt.df2(b,a);
hoverall=cascade(hallpass,hellip);
hFVT=fvtool([hellip,hoverall]);
set(hFVT,'Filter',[hellip,hoverall]);
legend(hFVT,'Lowpass Elliptic filter','Compensated filter');

clear;
[num1,den1]=ellip(5,0.8,35,0.35);
[GdH,w]=grpdelay(num1,den1,512);
plot(w/pi,GdH); grid
xlabel('\omega/\pi'); ylabel('Group delay, samples');
F=0:0.001:0.4;
g=grpdelay(num1,den1,F,2); % Equalize the passband
Gd=max(g)-g;

% Design the allpass delay equalizer
[num2,den2]=iirgrpdelay(10,F,[0,0.4],Gd);
[GdA,w] = grpdelay(num2,den2,512);
hold on;
plot(w/pi,GdH+GdA,'r');
legend('Original Filter','Compensated filter');
```

13. 编写4点滑动平均滤波器程序。原始未受干扰的序列为：s[n]=3n(0.8)^n，加性噪声信号d[n]为随机序列，幅度0.6，受干扰的序列为：x[n]= s[n]+ d[n]，分别绘制长度为40的原始未受干扰的序列，噪声序列和受干扰序列，以及滑动平均滤波器的输出。

```matlab
% Program 2_4
% Signal Smoothing by a Moving-Average Filter
R = 40;
d = 6/5*(rand(1,R)-0.5);
m = 0:1:R-1;
s =3.*m.*0.8.^m;
x = s + d;
subplot(211);
plot(m,d,'r-',m,s,'b:',m,x,'m--')
title('The sequence with noise');
ylabel('Amplitude')
legend('d[n]','s[n]','x[n]');
b = ones(4,1)/4;
y = fftfilt(b,x);
subplot(212);
plot(m,s,'r-',m,y,'b--')
title('The original sequence & the output sequence');
legend('s[n]','y[n]');
ylabel('Amplitude')
```

14. 绘制长度为10点的矩形序列的16点离散傅立叶变换样本的幅度和相位。

```matlab
xn=ones(10,1);
Xk=fft(xn,16);Xkf=abs(Xk);Xkp=angle(Xk);
subplot(211);
stem(0:15,Xkf,'filled');
xlabel('Time index/n');ylabel('Magnitude');
subplot(212);
stem(0:15,Xkp,'filled');
xlabel('Time index/n');ylabel('Phase')
```

15. 已知系统的系统函数为：

$$H(z) = \frac{z^2 - 0.2z + 0.5}{z^4 + 3.2z^3 + 1.5z^2 - 0.8z + 1.4}$$

用MATLAB的filter函数求h[n]的前20个样本值。

```matlab
num=[1,-0.2,0.5];
den=[1,3.2,1.5,-0.8,1.4];
x=[1 zeros(1,20-1)];
y=filter(num, den, x);
disp('Coefficients of the power series expansion');
disp(y)
stem(y)
```

16. 利用Hermann公式估计FIR低通滤波器的阶数。该滤波器的性能指标为：通带截止频率为1500Hz，阻带截止频率为1800Hz，通带波纹为0.015，阻带波纹为0.021，抽样频率为5000Hz。

```matlab
% Program 10_1
% Estimation of FIR Filter Order Using remezord
%
clear;
fedge = [1500 1800] ;% input('Type in the bandedges = ');
mval = [1 0];%input('Desired magnitude values in each band = ');
dev = [0.015 0.021];%input('Allowable deviation in each band = ');
FT = 5000;%input('Type in the sampling frequency = ');
[N, fpts, mag, wt] = remezord(fedge, mval, dev, FT);
d = fdesign.lowpass('n,fp,fst',N,0.6,0.72);
design(d);
fprintf('Filter order is %d \n',N);
```

17. 编写长度为5的中值滤波器程序。原始未受干扰的序列为：s[n]=3n(0.8)^n，加性噪声信号d[n]为随机序列，幅度0.6，分别绘制长度为40的受干扰序列，以及中值滤波器的输出。

```matlab
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
```

18. 已知16点序列x[n]的DFT为：

$$X[k] = \frac{k}{16}, k=0,1,...,15$$

绘制序列x[n]的实部和虚部。

```matlab
Xk=(0:15)/16;xn=ifft(Xk);
xnre=real(xn);
xnim=imag(xn);
subplot(2,1,1);
stem(0:15,xnre,'.');grid on;
title('The real part of the sequence');
subplot(2,1,2);
stem(0:15,xnim,'.');grid on;
title('The imaginary part of the sequence');
```

19. 已知系统的系统函数为：

$$H(z) = \frac{z^2 - 0.2z + 0.5}{z^4 + 3.2z^3 + 1.5z^2 + 0.8z + 1.4}$$

用MATLAB测试该系统的稳定性。

```matlab
num=[1 -0.2 0.5];
den=[1 3.2 1.5 0.8 1.4];
zplane(num,den);
grid on;
```

20. 利用Remez算法设计一个等波纹线性相位FIR低通滤波器。该滤波器的性能指标为：通带截止频率为1500Hz，阻带截止频率为1800Hz，通带波纹为0.015，阻带波纹为0.021，抽样频率为5000Hz。

```matlab
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
```

21. 已知序列$x_1[n]=[2.2, 3, -1.5, 4.2, -1.8]$，$x_2[n]=[0.8, -1, 1.6, 0.8]$，计算两个序列的卷积$y[n]=x_1[n]*x_2[n]$，并绘制序列$y[n]$。

```matlab
clc;
clear;
x1=[2.2 3 -1.5 4.2 -1.8];
x2=[0.8 -1 1.6 0.8];
x=conv(x1,x2) %结果在主界面输出
stem(x,'filled');
grid on;
xlabel('Time index/n');ylabel('Amplitude');
title('The output convolution');
```

22. 已知序列$x[n]=\cos(0.5\pi n)$，绘制序列x[n]的DFT和DTFT的幅度。

```matlab
n=0:15;x=cos(pi*n*0.5);
X=fft(x);
subplot(2,1,1);
stem(n,X,'.');
title('Magnitude of DFT')
xlabel('n'); ylabel('Magnitude')

%circulating DTFT
k=0:499;
w = pi/500*k;
X1=x*(exp(-1i*pi/500)).^(n'*k);
magX=abs(X1);
subplot(2,1,2);
plot(w/pi,magX);title('幅度响应');grid;
ylabel('幅度');xlabel('以\pi为单位的频率');
```

23. 已知FIR滤波器的系统函数为：

$$H(z) = 2.4 + 3.2z^{-1} + 1.5z^{-2} + 0.8z^{-3} + 1.4z^{-4} + 3.6z^{-5} + 5.2z^{-6}$$

用MATLAB将系统函数分解为二次多项式之积，并写出各二次多项式的表达式。

```matlab
clear
P=[2.4,3.2,1.5,0.8,1.4,3.6,5.2];
r=roots(P);%调用函数计算
syms z
s1=simple((z-r(1))*(z-r(2)));
d1=simple(s1./z^2)
s2=simple((z-r(3))*(z-r(4)));
d2=simple(s2./z^2)
s3=simple((z-r(5))*(z-r(6)));
d3=simple(s3./z^2)
Q=2.4*d1*d2*d3
```

24. 已知FIR数字低通滤波器的性能指标为：通带截止频率0.35π，阻带截止频率0.45π，通带和阻带波纹δ = 0.01。设计满足该滤波器的Kaiser's窗函数，绘制出Kaiser's窗函数的增益响应。

```matlab
clear;
fpts=[0.35,0.45];
mag=[1,0];
dev=[0.01,0.01];
[N,Wn,beta,ftype]=kaiserord(fpts,mag,dev);
kw=kaiser(N+1,beta);
b=fir1(N,Wn, kw);
[h,omega]=freqz(b,1,512);
plot(omega/pi,20*log10(abs(h)));grid;
xlabel('\omega/\pi'); ylabel('Gain, dB');
```

25. 已知系统的频响特性为：

$$H(e^{j\omega}) = \frac{1 - 0.2e^{-j\omega} + 0.5e^{-j2\omega} + 2e^{-j3\omega} - 0.6e^{-j4\omega}}{1 + 3.2e^{-j\omega} + 1.5e^{-j2\omega} - 0.8e^{-j3\omega} + 1.4e^{-j4\omega}}$$

绘制该系统的幅频特性和相频特性。

```matlab
clear
k=500; % number of frequency samples is 500
num=[1 -.2 .5 2 -.6]; %Numerator coefficients
den=[1 3.2 1.5 -.8 1.4]; %Denominator coefficients
w=0:pi/(k-1):pi;
h=freqz(num,den,w); %Compute the frequency response
subplot(1,2,1)
plot(w/pi,abs(h))
title('Magnitude Spectrum')
xlabel('\omega/\pi');ylabel('Magnitude')
subplot(1,2,2)
plot(w/pi,unwrap(angle(h))) %unwrapped phase function
title('Phase Spectrum')
xlabel('\omega/\pi');ylabel('Phase,radians')

[sos,g]=tf2sos(num,den)
```

26. 已知序列$x_1[n]=[2.2, 3, -1.5, 4.2, -1.8]$，$x_2[n]=[0.8, -1, 1.6, 0.8]$，基于DFT计算两个序列的卷积$y[n]=x_1[n]*x_2[n]$，并绘制基于DFT计算得到的$y[n]$。

```matlab
clc;
clear;
x1=[2.2 3 -1.5 4.2 -1.8];
x2=[0.8 -1 1.6 0.8];
n=size(x1,2);
m=size(x2,2);
X1=fft(x1,n+m-1);
X2=fft(x2,n+m-1);
X=X1.*X2;
x=ifft(X)
stem(x,'.');
```

27. 已知IIR滤波器的系统函数为：

$$H(z) = \frac{2 + 5z^{-1} + z^{-2} - 3z^{-3} + 4z^{-4} + 6z^{-5}}{1 + 3z^{-1} - 5z^{-2} + 2z^{-3} - 4z^{-4} + 3z^{-5}}$$

用MATLAB将系统函数表示为级联型结构形式，并写出各级联子系统的表达式。

```matlab
format short;
num = [2,5,1,-3,4,6];
den = [1,3,-5,2,-4,3];
[z,p,k]=tf2zp(num,den);
sos=zp2sos(z,p,k)
```

28. 用Kaiser's窗函数设计FIR数字高通滤波器，其滤波器的性能指标为：通带截止频率0.55π，阻带截止频率0.45π，通带和阻带波纹δ =0.04。绘制出该滤波器的增益响应。

```matlab
% Program 10_5
% Lowpass Filter Design Using the Kaiser Window
clc
%delta_s=50; %阻带衰减
delta_p=0.04; %通带衰减
%alpha_s=10^(-delta_s/20);
alpha_p=1-10^(-delta_p/20);
fpts = [0.45 0.55]; %截止频率1，截止频率2
mag = [0 1]; %截止频率1对应的幅度 截止频率2对应的幅度
dev = [alpha_p alpha_p]; %通带衰减 阻带衰减（线性）
[N,Wn,beta,ftype]=kaiserord(fpts,mag,dev);
b = fir1(N,Wn,'high',kaiser(N+1,beta));
[h,omega] = freqz(b,1,512);
plot(omega/pi,20*log10(abs(h)));grid;
xlabel('\omega/\pi'); ylabel('Gain, dB');
```

29. 绘制6点滑动平均滤波器的幅频特性和相频特性。

```matlab
h = ones(1,6)/6;
[H,w] = freqz(h, 1, 256);
subplot(211)
plot(w/pi, abs(H));
ylabel('Magnitude'); xlabel('\omega/\pi');
subplot(212)
ph = angle(H)*180/pi;
plot(w/pi,angle(H)*180/pi);
ylabel('Phase, degrees');xlabel('\omega/\pi');
```

30. 原始序列为：s[n]=3n(0.8)^n，加性噪声d[n]为随机序列，幅度0.6，受干扰的序列为：x[n]= s[n]+ d[n]，使用重叠相加法实现5点滑动平均滤波器对x[n]的处理。绘制未受干扰的序列s[n]和滤波器输出的有噪序列(利用fftfilt函数)。

```matlab
% Program 2_4
% Signal Smoothing by a Moving-Average Filter
R = 40;
d = 6/5*(rand(1,R)-0.5);
m = 0:1:R-1;
s =3.*m.*0.8.^m;
x = s + d;
subplot(211);
plot(m,d,'r-',m,s,'b:',m,x,'m--')
title('The sequence with noise');
ylabel('Amplitude')
legend('d[n]','s[n]','x[n]');
b = ones(5,1)/5;
y = fftfilt(b,x);
subplot(212);
plot(m,s,'r-',m,y,'b--')
title('The original sequence & the output sequence');
legend('s[n]','y[n]');
ylabel('Amplitude')
```

31. 已知IIR滤波器的系统函数为：

$$H(z) = \frac{2 + 5z^{-1} + z^{-2} - 3z^{-3} + 4z^{-4} + 6z^{-5}}{1 + 3z^{-1} - 5z^{-2} + 2z^{-3} - 4z^{-4} + 3z^{-5}}$$

用MATLAB对系统进行并联结构I型和并联结构II型分解。

```matlab
num = [2 5 1 -3 4 6];
den = [1 3 -5 2 -4 3];
[r1,p1,k1] = residuez(num, den); %系统并联I 型分解
[r2,p2,k2] = residue(num, den); %系统并联II 型分解
disp('Parallel Form I')
disp('Residues are'); disp(r1);
disp('Poles are at'); disp(p1);
disp('Constant value'); disp(k1);
disp('Parallel Form II')
disp('Residues are'); disp(r2);
disp('Poles are at'); disp(p2);
disp('Constant value'); disp(k2);

%结果出来后，寻找其中的共轭复数对。制成rri=[ ] ppi=[ ] 按照题目中系数，构造如下
rri=[r1(3), r1(4)];
ppi=[p1(3), p1(4)];
[b1,a1]=residuez(rri, ppi, 0);
disp('Parallel Form I 复数极点顶展开式分子系数'); disp(b1);
disp('Parallel Form I 复数极点顶展开式分母系数'); disp(a1);

rr2=[r2(3), r2(4)];
pp2=[p2(3), p2(4)];
[b2,a2]=residue(rr2, pp2, 0);
disp('Parallel Form II 复数极点顶展开式分子系数'); disp(b2);
disp('Parallel Form II 复数极点顶展开式分母系数'); disp(a2);
```

32. 用海明窗设计多频带FIR滤波器，该滤波器满足如下条件。在频率范围0到0.32π内幅度为0.6，在频率范围0.35π到0.65π内幅度为0.2，在频率范围0.68π到π内幅度为0.8。绘制出该滤波器的幅频特性。

```matlab
fpts = [0 0.32 0.35 0.65 0.68 1];
mval = [0.6 0.6 0.2 0.2 0.8 0.8];
b = fir2(100,fpts,mval);
[h,omega] = freqz(b,1,512);
plot(omega/pi,abs(h));grid;
xlabel('\omega/\pi'); ylabel('Magnitude');
```

34. 已知系统的系统函数为：

$$H(z) = \frac{z^4 - 0.2z^3 + 0.5z^2 + 2z - 0.6}{z^4 + 3.2z^3 + 1.5z^2 - 0.8z + 1.4}$$

绘制该系统的零极点分布图。

```matlab
方法一：num=[1 -0.2 0.5 2 -0.6];
den=[1 3.2 1.5 -0.8 1.4];
zplane(num,den)
```

```matlab
%34
%绘制该系统的零极点分布图
num=[1, -0.2, 0.5, 2, -0.6]; %分子系数降幂排列
den=[1, 3.2, 1.5, -0.8, 1.4]; %分母系数降幂排列
[z, p, k]=tf2zp(num, den);
zplane(z, p); %画零极点的函数，表示毫无压力
grid on;
title('系统零极点分布图');
```

35. 已知全通系统的系统函数为：

$$H(z) = \frac{1 + 0.4z^{-1} + 0.18z^{-2} - 0.2z^{-3}}{1 + 0.2z^{-1} - 0.18z^{-2} - 0.4z^{-3}}$$

用MATLAB求全通系统进行级联格型结构的乘法器系数。

```matlab
num=[1 0.4 0.18 -0.2];
>> k=poly2rc(num)
```

36. 已知有限长序列为：$x[n] = \sin(25\pi n/64)$，求该序列的64点离散傅立叶变换X[k]，绘制出X[k]的幅度。

```matlab
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
```

37. 设计4阶巴特沃兹模拟低通滤波器，其3-dB截止频率为1，绘制滤波器的增益响应。

```matlab
N = 4;
Wn = 1;
[num,den] = butter(N,Wn,'s');
[h,w] = freqs(num,den);
plot (w,20*log10(abs(h)));
xlabel('Frequency, Hz'); ylabel('Gain, dB');
title('The 4th-order IIR Butterworth Lowpass Filter ')
grid on
```

38. 已知系统的零极点分别如下：

零点：2.2, -1+i, -1-i, 1.4
极点：3.7+2i, 3.7-2i, -2.1-i, -2.1+i

求系统的系统函数H(z)。

```matlab
format long
zr =[2.2 -1+i -1-i 1.4];
pr =[3.7+2*i 3.7-2*i -2.1-i -2.1+i ];
% Transpose zero and pole row vectors
z = zr'; p = pr';
k = 1;
[num, den] = zp2tf(z, p, k);
disp('Numerator polynomial coefficients'); disp(num);
disp('Denominator polynomial coefficients'); disp(den);
```

39. 设计椭圆IIR数字低通滤波器，其性能指标为：通带截止频率为1000Hz，阻带截止频率为1250Hz，通带波纹为0.4dB，最小阻带衰减为45dB，抽样频率为5000Hz。绘制所设计的滤波器增益响应。

```matlab
Fp = 1000%input('passband edge in Khz = ');
Fs = 1250%input('stopband edge in Khz = ');
Ft = 5000%input('Sampling rate in Khz = ');
Rp =0.4% input('Passband ripple in dB = ');
Rs =45% input('Minimum stopband attenuation in dB = ');
Wp=2*Fp/Ft;
Ws=2*Fs/Ft;
[N,Wn] = ellipord(Wp,Ws,Rp,Rs);
[b,a] = ellip(N,Rp,Rs,Wn);
[h,omega] = freqz(b,a,256);
plot (omega/pi,20*log10(abs(h)));grid;
xlabel('\omega/\pi'); ylabel('Gain, dB');
title('IIR Elliptic Lowpass Filter');
%figure(2);
%subplot(2,1,1);
%plot(omega/pi,20*log10(abs(h))); grid;
%axis([0 1 -60 5]);
%subplot(2,1,2);
%plot(omega/pi,20*log10(abs(h))); grid;
%axis([0 0.4 -0.6 0.2]);
```

40. 编写总体均值滤波器程序。原始未受干扰的序列为：s[n]=3n(0.8)^n，加性噪声信号d[n]为随机序列，幅度0.6，受干扰的序列为：x[n]= s[n]+ d[n]，绘制噪声序列和60次检测结果的总体平均的序列。

```matlab
% Program 2_4
% Signal Smoothing by a Moving-Average Filter
R = 60;
d = 6/5*(rand(1,R)-0.5);
m = 0:1:R-1;
s =3.*m.*0.8.^m;
x = s + d;
subplot(211);
plot(m,d,'r-',m,s,'b:',m,x,'m--')
title('The sequence with noise');
ylabel('Amplitude')
legend('d[n]','s[n]','x[n]');
b = ones(R,1)/R;
y = fftfilt(b,x);
subplot(212);
plot(m,s,'r-',m,y,'b--')
title('The original sequence & the output sequence');
```