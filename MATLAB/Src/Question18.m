Xk=(0:15)/16;xn=ifft(Xk);
xnre=real(xn);
xnim=imag(xn);
subplot(2,1,1);
stem(0:15,xnre,'.');grid on;
title('The real part of the sequence');
subplot(2,1,2);
stem(0:15,xnim,'.');grid on;
title('The imaginary part of the sequence');
