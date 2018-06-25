function [I sol_data] = runModelCholera(x, time_mesh, y0, solver)
% Solves the cholera model, with parameter perturbation indicated by 
% random vector x, and returns the solution values at time vector time_mesh

%
% Set the output flag (if set to 1, will show a plot of solution, 
%    and show some solver statistics).
%
output_flag = 0;

mu = [1.5;            
      7.5;
      1e6;
      1 / 1560; 
      168 / 5; 
      70;
      7/30; 
      7/5];

% Note: mu is the nominal paramater vector. 
% Elements of mu are as follows (coefficients
% are as in paper of Hartley et al.):
%    par(1) = beta_L
%    par(2) = beta_H
%    par(3) = kappa_L
%       Note: kappa_H = par(3)/700
%    par(4) = b 
%    par(5) = chi 
%    par(6) = xi 
%    par(7) = delta_L
%    par(8) = gamma

alpha = 0.1;
params = mu + alpha * mu .* x(:); 


if nargin == 3
   solver = @ode45; 
   N = y0(1)+y0(2)+y0(3);
elseif nargin == 2
   N = 1e4;
   y0 = [N-1 1 0 0 0];    % default initial condition
   solver = @ode45;
elseif nargin < 2
   error('Missing input arguments.');
end

% call solver
[sol_data] = modelSolveCholera(params, time_mesh, y0, solver);
I = sol_data.I;

%
% output
% 
if output_flag == 1
   stats = sol_data.stats;
   disp('Solver performance statistics:');
   disp(stats);
   close all; 
   loglog(time_mesh, I, 'linewidth', 2);
   xlabel('time (weeks)');
   ylabel('infected population');
   set(gca, 'fontsize', 20);
   xlim([1e-2 time_mesh(end)]);

   %hold on;
   %ts = sol_data.ts;
   %Is = sol_data.Is;
   %loglog(ts, Is, 'o', 'linewidth', 2);
   %xlabel('time (weeks)');
   %ylabel('infected population');
   %set(gca, 'fontsize', 20);
   %xlim([1e-2 time_mesh(end)]);
end




%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% Subfunction

function [sol_data] = modelSolveCholera(par, time_mesh, y0, ode_solver)
global N beta_L beta_H kappa_L kappa_H b chi xi delta_L gamma_I

% Note: I rename the parameter gamma to gamma_I for clarity and
%       also to avoid conflict with Matlab function gamma
beta_L = par(1); 
beta_H = par(2);
kappa_L = par(3); 
kappa_H = par(3) / 700; 
b = par(4); 
chi = par(5); 
xi = par(6);
delta_L = par(7); 
gamma_I = par(8); 
N = 10000;

%fprintf('%10.6e\n', [beta_L beta_H kappa_L kappa_H b chi xi delta_L gamma_I]);

%options = odeset('NonNegative', [1:5], 'abstol', 1e-6, 'reltol', 1e-6);
options = odeset('abstol', 1e-6, 'reltol', 1e-6);
Tf = time_mesh(end);
sol = ode_solver(@frhs, [0 Tf], y0, options); 

% Note: For MC runs, can comment out the S, R, BH, BL lines
%sol_data.S = deval(sol, time_mesh, 1);
sol_data.I = deval(sol, time_mesh, 2);
%sol_data.R = deval(sol, time_mesh, 3);
%sol_data.BH = deval(sol, time_mesh, 5);
%sol_data.BL = deval(sol, time_mesh, 4);

sol_data.Is = sol.y(2,:)';
sol_data.ts = sol.x;
sol_data.stats = sol.stats;



function dydt = frhs(t, y)
global N beta_L beta_H kappa_L kappa_H b chi xi delta_L gamma_I

S = y(1);
I = y(2);
R = y(3);
BH = y(4);
BL = y(5);

FL = BL / (kappa_L + BL);
FH = BH / (kappa_H + BH);
dydt = [         b * N - beta_L * S * FL - beta_H * S * FH - b * S;
                 beta_L * S * FL + beta_H * S * FH - (gamma_I + b) * I;
                 gamma_I * I - b * R;
                 xi * I - chi * BH;
                 chi * BH - delta_L * BL];
