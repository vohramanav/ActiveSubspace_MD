function Fhat = get_polyfit_surrogate(y, f, eta, deg)

if nargin < 3
   deg = 2;
end

p = polyfit(y(:), f(:), deg);

Fhat = @(x)(polyval(p,eta'*x));
