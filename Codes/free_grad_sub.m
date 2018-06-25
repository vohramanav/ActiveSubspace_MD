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

N = nsams - 1;
p = N - 1;
d_matrix = zeros(N,1);
M = 30;
% model realizations
G = zeros(nsams-1,1); 
G(1:23,1) = k_data(1:23,2);
G(24:nsams-1,1) = k_data(25:nsams,2);

% Draw M independent samples
samples = zeros(M,ncols); ypr = zeros(M,ncols);
for j = 1:dim
  samples(:,j) = unifrnd(L(j,1),U(j,1),M,1);
end

% project samples in [-1,1]
for i = 1:M
  for j = 1:dim
    ypr(i,j) = 2.0.*(samples(i,j)-L(j))./(U(j)-L(j)) - 1.0;
  end
end

for i = 1:M
  for j = 1:N
    d_matrix(j) = 0;
    for k = 1:dim
      d_matrix(j) = d_matrix(j) + norm(ypr(i,k) - xpr(j,k));
    end
  end
  [z,index(i,:)] = sort(d_matrix);
  for j = 1:p
    ip = (i-1)*p + j;
    points(ip,:) = xpr(index(i,j),:);
  end
end

% formulate the least squares problem

for np = 1:M
  A = [1 points((np-1)*p+1,:)];
  
  for i = (np-1)*p+2 : np*p        
    A = [A; 1 points(i,:)];    
  end

  B = G(index(np,1));

  for i=2:p
    B = [B; G(index(np,i))];
  end

  z = A \ B;
 
  if np == 1
    b_matrix = z(2:dim+1); %remove first entry of z 
  else
    b_matrix = [b_matrix z(2:dim+1)];
  end

end

% construct the covariance matrix

C = 0;

for i=1 : M
  z = b_matrix(:,i);
  C = C + z * z';
end

C=C/M;
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

eigv(lambda);
plot_gsa(as);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%[rn1,xmesh_PCE1,density_PCE1,gsa_m1,gsa_t1] = init_1D(eta,xpr,G,dim,L,U); % 1D 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
print -depsc free_eigv.eps
end

function gsa = plot_gsa(as);
gsa_t1 = load('free_gsa_t_1D.txt');
TI = zeros(7,2);
TI(:,1) = as;
TI(:,2) = gsa_t1;
figure;
bar(TI,'BarWidth',0.5);
%ylabel('$$\mathrm{Activity~Scores}$$','interpreter','latex','fontsize',18);
xtickl = ({'$$\mathrm{A}$$','$$\mathrm{B}$$','$$\mathrm{p}$$','$$\mathrm{q}$$',...
           '$$\mathrm{\alpha}$$','$$\mathrm{\lambda}$$','$$\mathrm{\gamma}$$',...
           'interpreter','latex'});
leg = legend({'$\mathrm{\mathcal{T}(\theta_i)}$','$\mathrm{\eta_i}$'});
set(leg,'interpreter','latex','fontsize',16,'location','NorthWest');
set(gca,'xtick',1:7,'xticklabel',xtickl,'fontsize',16);
set(gca,'TickLabelInterpreter','latex');
set(gcf,'color',[1,1,1]);
print -depsc free_as_gsa.eps

end

























