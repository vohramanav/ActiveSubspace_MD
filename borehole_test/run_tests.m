close all;
clear all;

% gradient based
active_subspace;


% gradient free
local_linear_approx;


% compare dominant eigenvectors
figure;
hold on;

% line them up: 
V(:,1) = V(:,1)*sign(V(1,1));
W(:,1) = W(:,1)*sign(W(1,1));

plot(V(:,1), '-o', 'linewidth',2);
plot(W(:,1), '-*', 'linewidth',2);

legend('alg 1.1', 'alg 1.2');
set(gca, 'fontsize', 20);
title('comparing dominant eigenvectors');


% compare normalized eigenvalues
l1 = lambda_grad./lambda_grad(1);
l2 = lambda_loclin./lambda_loclin(1);
figure;
bar([l1(:) l2(:)]);
title('comparing eigenvalues');
set(gca, 'yscale', 'log');
set(gca, 'fontsize', 20);
legend('alg 1.1', 'alg 1.2');

