clear
P = [2.4, 3.2, 1.5, 0.8, 1.4, 3.6, 5.2];
r = roots(P);
syms z

% 直接计算，不使用 simple
d1 = expand((z - r(1)) * (z - r(2))) / z^2;
d2 = expand((z - r(3)) * (z - r(4))) / z^2;
d3 = expand((z - r(5)) * (z - r(6))) / z^2;

Q = 2.4 * d1 * d2 * d3;
pretty(Q)  % 以更美观的形式显示