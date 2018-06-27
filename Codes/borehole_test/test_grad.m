ndim = 7;
% base point
xi = -1 + 2 * rand(ndim, 1);

[f0 df] = borehole(xi);

% compute FD directional derivative

% random direction
d = randn(ndim, 1); d = d / norm(d);

% analytic directional derivative
deriv_analytic = dot(df, d);

t = 1;

fprintf('  t               rel_err        ratio\n');
fprintf('--------------------------------------\n');
for k = 1 : 10
   t = t / 2;

   f1 = borehole(xi + t * d);
   deriv_fd(k) = (f1 - f0) / t;

   err(k) = abs(deriv_fd(k) - deriv_analytic) / abs(deriv_analytic);
   if k == 1
      fprintf('%10.6f  %16.4e     .....\n', t, err(k));
   else
      ratio = err(k-1)/err(k);
      fprintf('%10.6f  %16.4e %10.4f\n', t, err(k), ratio);
   end
end

