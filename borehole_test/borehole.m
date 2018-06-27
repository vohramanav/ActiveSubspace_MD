function [f Df] = borehole(xi)
%
%   Input: 
%      xi --- 7x1 vector with entries in [-1, 1]
%
%   Output:
%
%      f --- function values
%     Df --- 7x1 gradient vector


% define the physical ranges
L = [0.05 1120 63070 990 63 700 9855]';
U = [0.15 1680 115600 1110 1116 820 12045]';
 
% map x to physical ranges
x = 0.5 * (L+U) + 0.5 * (U-L) .* xi(:); 

x1 = x(1);
x2 = x(2);
x3 = x(3);
x4 = x(4);
x5 = x(5);
x6 = x(6);
x7 = x(7);

r = 3698.30;

% function eval
f = 2*pi*x3*(x4-x6) / (log(r/x1)*(1+2*x2*x3/(log(r/x1)*x1^2*x7)+x3/x5));

% derivative
Df_dx(1) = (2*(x4-x6))*x3*x1*pi*x7*x5*(x7*(x3+x5)*x1^2+4*x2*x3*x5)/(x7*(x3+x5)*x1^2*log(r/x1)+2*x2*x3*x5)^2;
Df_dx(2) = -4*pi*x3^2*(x4-x6)*x1^2*x7*x5^2/(x7*(x3+x5)*x1^2*log(r/x1)+2*x2*x3*x5)^2;
Df_dx(3) = 2*pi*(x4-x6)*x1^4*x7^2*x5^2*log(r/x1)/(x7*(x3+x5)*x1^2*log(r/x1)+2*x2*x3*x5)^2;
Df_dx(4) = 2*pi*x3*x1^2*x7*x5/(x7*(x3+x5)*x1^2*log(r/x1)+2*x2*x3*x5);
Df_dx(5) = 2*pi*x3^2*(x4-x6)*log(r/x1)*x1^4*x7^2/(x7*(x3+x5)*x1^2*log(r/x1)+2*x2*x3*x5)^2;
Df_dx(6) = -2*pi*x3*x1^2*x7*x5/(x7*(x3+x5)*x1^2*log(r/x1)+2*x2*x3*x5);
Df_dx(7) = (4*pi*x3^2*(x4-x6)*x1^2*x5^2*x2)/((x7*(x3+x5)*x1^2*log(r/x1)+2*x2*x3*x5)^2);

% chain rule
dx_dxi = 0.5 * (U-L);

Df = Df_dx(:) .* dx_dxi; 
