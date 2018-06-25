close all;
clear all;
rng(100);
set(0,'DefaultFigureVisible','off');

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

nsams = 35;
% refine xp as xpr
xpr = zeros(nsams-1,ncols);
xpr(1:23,:) = xp(1:23,:);
xpr(24:nsams-1,:) = xp(25:nsams,:);

for j = 1:dim
  dX(j) = 1e-5.*(U(j,1)-L(j,1));
end

G = zeros(nsams-1,1); Gdx = zeros(nsams-1,dim);
G(1:23,1) = k_data(1:23,2);
G(24:nsams-1,1) = k_data(25:nsams,2);
dGdx = zeros(nsams-1,dim);

for j = 1:dim
  Gdx(1:23,j) = k_data(j*nsams+1:j*nsams+23,2);
  Gdx(24:nsams-1,j) = k_data(j*nsams+25:j*nsams+nsams,2);
  dGdx(:,j) = (Gdx(:,j) - G(:,1))./dX(j);
end

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


% Computing activity scores

as = zeros(7,1);

for i = 1:dim
  for j=1:3
    as(i) = as(i) + lambda(j).*(W(i,j).^2);
  end
end

as = as./sum(as);

% Eigvalue plot
%eigv(lambda);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%[rn1,xmesh_PCE1,density_PCE1,gsa_m1,gsa_t1] = init_1D(eta,xpr,G,dim,L,U); % 1D active var analysis, 
%[rn2,xmesh_PCE2,density_PCE2,gsa_m2,gsa_t2] = init_2D(eta,xpr,G,dim,L,U,xmesh_PCE1,density_PCE1); 
% 2D active var analysis, rn: rel L2 norm of error
%rn3 = init_3D(eta,xpr,G,dim,L,U,xmesh_PCE1,density_PCE1,xmesh_PCE2,density_PCE2); % 2D active var analysis
%[rn3,gsa_m3,gsa_t3] = init_3D(eta,xpr,G,dim,L,U); % 2D active var analysis

plot_gsa(as);
%plot_rn(rn1,rn2,rn3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function eigen = eigv(lambda)
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
print -depsc eigv.eps
end

function prn = plot_rn(rn1,rn2,rn3);
figure;
plot([1 2 3],[rn1 rn2 rn3],'o--','MarkerFaceColor','k');
xlabel('$$\mathrm{Active~Subspace~Dimension}$$','interpreter','latex','fontsize',18);
ylabel('$$\mathrm{Relative~L-2~Norm~of~Error}$$','interpreter','latex','fontsize',18);
xticks([1 2 3]);
xlim([0,4]);
set(gca,'TickLabelInterpreter','latex');
set(gcf,'color',[1,1,1]);
box on;
print -depsc rn.eps
end

function gsa = plot_gsa(as);
gsa_t1 = load('gsa_t_1D.txt');
gsa_t2 = load('gsa_t_2D.txt');
gsa_t3 = load('gsa_t_3D.txt');
TI = zeros(7,3);
TI(:,1) = gsa_t1;
TI(:,2) = gsa_t2;
TI(:,3) = gsa_t3;

figure;
bar(TI);
ylabel('$$\mathrm{Sobol~Total~Effect~Index}$$','interpreter','latex','fontsize',18);
xtickl = ({'$$\mathrm{A}$$','$$\mathrm{B}$$','$$\mathrm{p}$$','$$\mathrm{q}$$',...
           '$$\mathrm{\alpha}$$','$$\mathrm{\lambda}$$','$$\mathrm{\gamma}$$',...
           'interpreter','latex'});
leg = legend({'$\mathrm{1D~Subspace}$','$\mathrm{2D~Subspace}$','$\mathrm{3D~Subspace}$'});
set(gca,'xtick',1:7,'xticklabel',xtickl,'fontsize',16);
set(leg,'interpreter','latex','fontsize',14);
set(gca,'TickLabelInterpreter','latex');
set(gcf,'color',[1,1,1]);
print -depsc gsa_subspace.eps

figure;
bar(as,'BarWidth',0.5);
ylabel('$$\mathrm{Activity~Scores}$$','interpreter','latex','fontsize',18);
xtickl = ({'$$\mathrm{A}$$','$$\mathrm{B}$$','$$\mathrm{p}$$','$$\mathrm{q}$$',...
           '$$\mathrm{\alpha}$$','$$\mathrm{\lambda}$$','$$\mathrm{\gamma}$$',...
           'interpreter','latex'});
set(gca,'xtick',1:7,'xticklabel',xtickl,'fontsize',16);
set(gca,'TickLabelInterpreter','latex');
set(gcf,'color',[1,1,1]);
print -depsc activity_scores.eps

end

