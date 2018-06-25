function [rn1,xmesh_PCE1,density_PCE1,gsa_m,gsa_t] = init_1D(eta,xpr,G,dim,L,U);

% 1D Surrogate
y = eta(:,1)'*xpr'; % active var
Y_surr = get_polyfit_surr(y,G,eta(:,1),2);

% Compute rel L-2 norm of error at MD data for kappa
np = size(xpr,1); % number of points at which MD data is available
Y_surr_data = zeros(np,1);

for i = 1:np
  Y_surr_data(i,1) = Y_surr(xpr(i,:)');
end

Nr = (sum((G-Y_surr_data).^2)).^(0.5);
Dr = (sum((G).^2)).^(0.5);
rn1 = Nr./Dr;


% SSP 1D
ssp1D(eta,xpr,G);

% Generate a million random samples in the 7D space
ts = 1e5;
samples = zeros(ts,dim);

for j = 1:dim
  samples(:,j) = unifrnd(L(j,1),U(j,1),ts,1);
end

% project those samples in [-1,1]

xrand = zeros(ts,dim);

for i = 1:ts
  for j = 1:dim
    xrand(i,j) = 2.0.*(samples(i,j)-L(j))./(U(j)-L(j)) - 1.0;
  end
end

% Estimate k using 1D surrogate for the million projected samples
Y_surr_rs = zeros(ts,1);

for i = 1:ts
  Y_surr_rs(i,1) = Y_surr(xrand(i,:)');
end

% compare pdf of kappa using 1D surrogate
pdf_comp_ssp1D(G,Y_surr_rs);
[gsa_m,gsa_t] = sobol(L,U,dim);

% Function Definitions

% 1D polyfit
function surr = get_polyfit_surr(y,G,eta,deg)
p = polyfit(y(:),G,deg);
surr = @(x)(polyval(p,eta'*x));
end

function [gsa_m,gsa_t] = sobol(L,U,dim)
N = 1e5; % number of samples
A = zeros(N,dim); B = zeros(N,dim);

for i = 1:dim
  A(:,i) = unifrnd(L(i,1),U(i,1),N,1);
  B(:,i) = unifrnd(L(i,1),U(i,1),N,1);
end

% project those samples in [-1,1]

An = zeros(N,dim); Bn = zeros(N,dim);

for j = 1:N
  for i = 1:dim
    An(j,i) = 2.0.*(A(j,i)-L(i))./(U(i)-L(i)) - 1.0;
    Bn(j,i) = 2.0.*(B(j,i)-L(i))./(U(i)-L(i)) - 1.0;
  end
end

% constructing the modified matrices A(B)i and B(A)i
An_Bn_1 = An; An_Bn_1(:,1) = Bn(:,1);
An_Bn_2 = An; An_Bn_2(:,2) = Bn(:,2);
An_Bn_3 = An; An_Bn_3(:,3) = Bn(:,3);
An_Bn_4 = An; An_Bn_4(:,4) = Bn(:,4);
An_Bn_5 = An; An_Bn_5(:,5) = Bn(:,5);
An_Bn_6 = An; An_Bn_6(:,6) = Bn(:,6);
An_Bn_7 = An; An_Bn_7(:,7) = Bn(:,7);

Bn_An_1 = Bn; Bn_An_1(:,1) = An(:,1);
Bn_An_2 = Bn; Bn_An_2(:,2) = An(:,2);
Bn_An_3 = Bn; Bn_An_3(:,3) = An(:,3);
Bn_An_4 = Bn; Bn_An_4(:,4) = An(:,4);
Bn_An_5 = Bn; Bn_An_5(:,5) = An(:,5);
Bn_An_6 = Bn; Bn_An_6(:,6) = An(:,6);
Bn_An_7 = Bn; Bn_An_7(:,7) = An(:,7);

%surrogate evaluations
fAn = zeros(N,1); 
fAn_Bn = zeros(N,7); fBn_An = zeros(N,7);

for i = 1:N
  fAn(i,1) = Y_surr(An(i,:)');
  fAn_Bn(i,1) = Y_surr(An_Bn_1(i,:)'); fBn_An(i,1) = Y_surr(Bn_An_1(i,:)'); 
  fAn_Bn(i,2) = Y_surr(An_Bn_2(i,:)'); fBn_An(i,2) = Y_surr(Bn_An_2(i,:)'); 
  fAn_Bn(i,3) = Y_surr(An_Bn_3(i,:)'); fBn_An(i,3) = Y_surr(Bn_An_3(i,:)'); 
  fAn_Bn(i,4) = Y_surr(An_Bn_4(i,:)'); fBn_An(i,4) = Y_surr(Bn_An_4(i,:)'); 
  fAn_Bn(i,5) = Y_surr(An_Bn_5(i,:)'); fBn_An(i,5) = Y_surr(Bn_An_5(i,:)'); 
  fAn_Bn(i,6) = Y_surr(An_Bn_6(i,:)'); fBn_An(i,6) = Y_surr(Bn_An_6(i,:)'); 
  fAn_Bn(i,7) = Y_surr(An_Bn_7(i,:)'); fBn_An(i,7) = Y_surr(Bn_An_7(i,:)'); 
end

f_total = zeros((dim+1).*N,1);
f_total(1:N,1) = fAn;
f_total(N+1:end,1) = reshape(fBn_An,N*dim,1);
f0_m = mean(f_total); Dr_m = var(f_total);
f_total(N+1:end,1) = reshape(fAn_Bn,N*dim,1);
f0_t = mean(f_total); Dr_t = var(f_total);

Nr_m = zeros(dim,1); Nr_t = zeros(dim,1); gsa_m = zeros(dim,1); gsa_t = zeros(dim,1);

for i = 1:dim
   Nr_m(i) = (1.0./N).*dot(fAn,fBn_An(:,i)) - f0_m.^2.0;
%   Nr_t(i) = (1.0./N).*dot(fAn,fAn_Bn(:,i)) - f0_t.^2.0;
   diff = fAn - fAn_Bn(:,i);
   Nr_t(i) = (1.0/N).*dot(fAn,diff); 
end

gsa_m = Nr_m./Dr_m;
%gsa_t = 1.0 - Nr_t./Dr_t;
gsa_t = Nr_t./Dr_t;

end

% plot 1D SSP
function ssp1 = ssp1D(eta,xpr,G)
figure(1)
y = eta(:,1)'*xpr';
plot(y,G,'ko','markerfacecolor','k');
xlabel('$$\mathrm{\eta_1^{T}\xi}$$','interpreter','latex','fontsize',18);
ylabel('$$\mathrm{\kappa}$$','interpreter','latex','fontsize',18);
set(gca,'TickLabelInterpreter','latex','fontsize',14);
set(gcf,'color',[1,1,1]);
box on;
print -depsc free_ssp1D.eps
end

% compare pdf of kappa using 1D SSP
function comp1 = pdf_comp_ssp1D(G,Y_PCE)
step = (max(Y_PCE) - min(Y_PCE)).*(1e-4);
pts_PCE = min(Y_PCE):step:max(Y_PCE);
[density_PCE1,xmesh_PCE1] = ksdensity(Y_PCE,pts_PCE);

figure;
hold on;
nbins = 15;
histogram(G,nbins,'Normalization','pdf','FaceAlpha',0.2);
plot(xmesh_PCE1,density_PCE1,'Linewidth',2,'color','b');
%plot(xmesh_Model,density_Model,'Linewidth',2,'color','k');
xlabel('$$\mathrm{\kappa}$$','interpreter','latex');
ylabel('$$\mathrm{PDF}$$','interpreter','latex');
set(gca,'TickLabelInterpreter','latex','fontsize',14);
set(gcf,'color',[1,1,1]);
leg = legend({'$\mathrm{MD~Prediction}$','$\mathrm{1D~Subspace}$'});
set(leg,'Interpreter','latex');
box on;
print -depsc free_pdf_comp_SSP1D.eps
end

end
