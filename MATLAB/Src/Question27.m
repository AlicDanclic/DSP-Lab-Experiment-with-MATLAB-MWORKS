format short;
num = [2,5,1,-3,4,6];
den = [1,3,-5,2,-4,3];
[z,p,k]=tf2zp(num,den);
sos=zp2sos(z,p,k)
