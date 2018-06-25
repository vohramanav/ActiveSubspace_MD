close all;
clear all;
rng(100);

dim = 7;
dX = zeros(dim,1); 
L = zeros(dim,1); U = zeros(dim,1);
k_data = load('k_data35.txt');
params = load('params_35.txt');

% Nominal values of the parameters
A = 7.049556277; B = 0.6022245584;
p = 4.0; q = 0.0; alpha = 1.80;
lambda = 21.0; gamma = 1.20;

N = [A;B;p;q;alpha;lambda;gamma];
L(:,1) = 0.9.*N(:); % lower-bound
U(:,1) = 1.1.*N(:); % upper-bound
U(4,1) = 0.1;

% Project params in [-1,1]

nrows = size(params,1); % number of rows
ncols = size(params,2); % number of cols

xp = zeros(nrows,ncols);

for i = 1:nrows
  for j = 1:ncols
    xp(i,j) = 2.0.*(params(i,j)-L(j))./(U(j)-L(j)) - 1.0;
  end
end


for j = 1:dim
  dX(j) = 1e-5.*(U(j,1)-L(j,1));
end

nsams = 35;
G = zeros(nsams-1,1); Gdx = zeros(nsams-1,dim);
G(1:23,1) = k_data(1:23,2);
G(24:nsams-1,1) = k_data(25:nsams,2);
dGdx = zeros(nsams-1,dim);

for j = 1:dim
  Gdx(1:23,j) = k_data(j*nsams+1:j*nsams+23,2);
  Gdx(24:nsams-1,j) = k_data(j*nsams+25:j*nsams+nsams,2);
  dGdx(:,j) = (Gdx(:,j) - G(:,1))./dX(j);
end

% refine xp as xpr

xpr = zeros(nsams-1,ncols);
xpr(1:23,:) = xp(1:23,:);
xpr(24:nsams-1,:) = xp(25:nsams,:);

% Computing the C matrix (Eq. 1.2)

C = zeros(dim,dim);
grad_f = zeros(dim,1);

for i = 1:nsams-1
  grad_f(:,1) = dGdx(i,:);
  C = C + grad_f*transpose(grad_f);
end

C = (1./(nsams-1)).*(C);
[W,D] = eig(C);

[lambda, idx] = sort(diag(D), 'descend');
W = W(:,idx);

eta = zeros(dim,3);
eta(:,1) = W(:,1);
eta(:,2) = W(:,2);
eta(:,3) = W(:,3);

% Activity Scores
%ac = activity_scores(dim,lambda,eta);

% Surrogates
y = eta(:,1)'*xpr'; % active var
Y_surr = get_polyfit_surr(y,G,eta(:,1),3);

% Distribution for kappa using the polyfit surrogate

% Generate samples in the original 7D space
ts = 1e6;
samples = zeros(ts,dim);

for j = 1:dim
  samples(:,j) = unifrnd(L(j,1),U(j,1),ts,1);
end

% project those samples in [-1,1]

% compute the distribution using the 1D surrogate and plot it 


% Plots
%eig(lambda);
ssp1D(eta,xpr,G);
%ssp2D(eta,xpr,G);

function as = activity_scores(dim,lambda,eta)
ac = zeros(dim,1);

for i = 1:3
  ac(:) = ac(:) + lambda(:).*(eta(:,i).^2);
end

ac(:) = ac(:)./max(ac);
end

function surr = get_polyfit_surr(y,G,eta,deg)
p = polyfit(y(:),G,deg);
surr = @(x)(polyval(p,eta'*x));
end

% plot of eigenvalues
function eigen = eig(lambda)
figure(1)
%plot(1:7,e,'--o','MarkerFaceColor','b');
semilogy(lambda./lambda(1),'-o');
xlabel('$$\mathrm{index}$$','interpreter','latex','fontsize',18);
ylabel('$$\mathrm{Eigenvalues}$$','interpreter','latex','fontsize',18);
set(gca, 'xtick',1:7,'fontsize',14);
set(gca,'TickLabelInterpreter','latex');
set(gcf,'color',[1,1,1]);
box on;
grid on;
print -depsc eig.eps
end

%% univariate SSP
function ssp1 = ssp1D(eta,xpr,G)
figure(1)
y = eta(:,1)'*xpr';
plot(y,G,'ko','markerfacecolor','k');
xlabel('$$\mathrm{\eta^{T}\xi}$$','interpreter','latex','fontsize',18);
ylabel('$$\mathrm{\kappa}$$','interpreter','latex','fontsize',18);
set(gca,'TickLabelInterpreter','latex','fontsize',14);
set(gcf,'color',[1,1,1]);
box on;
print -depsc ssp1D.eps
end


%% bivariate SSP
function ssp2 = ssp2D(eta,xpr,G)
figure(2)
y1 = eta(:,1)'*xpr';
y2 = eta(:,2)'*xpr';
scatter3(y1,y2,G,'ko','markerfacecolor','k');
xlabel('$$\mathrm{\eta_1^{T}\xi}$$','interpreter','latex','fontsize',18);
ylabel('$$\mathrm{\eta_2^{T}\xi}$$','interpreter','latex','fontsize',18);
zlabel('$$\mathrm{\kappa}$$','interpreter','latex','fontsize',18);
set(gca,'TickLabelInterpreter','latex','fontsize',14);
set(gcf,'color',[1,1,1]);
box on;
print -depsc ssp2D.eps
end




