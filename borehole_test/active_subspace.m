m = 7;
k = m + 1;
alpha = 8;
rng(1245);
nsamples = floor(alpha * k * log(m));

x = -1 + 2 * rand(nsamples, m);

fprintf('running %i samples\n', nsamples);
parfor i = 1 : nsamples
   fprintf('%i\n', i);
   [y(i) dy(i,:)] = borehole(x(i,:)');
end

C = 0;

for i = 1 : nsamples
   dfi = dy(i,:)';
   C = C + dfi * dfi'; 
end
C = C / nsamples;

[V D] = eig(C);
[lambda_grad, idx] = sort(diag(D), 'descend');
V = V(:,idx);

figure
semilogy(abs(lambda_grad)./lambda_grad(1), '-o');
set(gca, 'fontsize', 20);
title('grad based eigs');


eta1 = V(:,1);
eta2 = V(:,2);

%
% univariate
%
figure;
g1 = eta1'*x';
plot(g1, y, 'ko', 'markerfacecolor', 'k')
set(gca, 'fontsize', 20);
xlabel('<eta1, x>');
ylabel('f(x)');
title('grad based SSP');







