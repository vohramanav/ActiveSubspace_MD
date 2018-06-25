F = @(x)(cholera_scalar_qoi(x));

ndim = 8;
dist = 'uniform';

% get the active subspace 
[lambda, W, X, ff, yy] = compute_active_subspace(F, ndim, dist);

% construct a surrogate using the active subspace
Fhat = get_polyfit_surrogate(yy, ff, W(:,1), 2);

Ns = length(X);

for i = 1 : Ns
   ff_hat(i) = Fhat(X(i,:)');
end

Ns = 1000;
X_rand = -1 + 2 * rand(Ns, ndim);

for i = 1 : Ns
   ff_hat_samp(i) = Fhat(X_rand(i,:)');
end

[pdf xi] = ksdensity(ff_hat_samp);

%
% figures
%
%close all;
figure(1);
plot(ff, ff_hat, '*');
set(gca, 'fontsize', 20);
xlabel('true model');
ylabel('surrogate');
print -depsc surr.eps

figure(2)
hold on;
histnorm(ff);
plot(xi, pdf, 'r', 'linewidth',2);
xlabel('parameter x');
ylabel('distribution');
legend('true model evals', 'surrogate');
set(gca, 'fontsize', 20);
print -depsc pdf_comp.eps
