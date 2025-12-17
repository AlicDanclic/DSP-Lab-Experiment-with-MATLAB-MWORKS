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
