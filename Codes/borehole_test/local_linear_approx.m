m = 7;
F = @(x)(borehole(x));

%Algorithm 1.2 from book

k = m + 1;
alpha = 8;

N = alpha * m + 1;
M = floor( alpha *k * log(m) );

%
% get the samples
%
x = -1 + 2 * rand(N , m); 
y = -1 + 2 * rand(M , m); 

fprintf('running %i samples\n', N);

%
% Model runs
%
parfor i = 1 : N
    
   fprintf('%i\n', i);
   
   [f(i),~] = F(x(i,:)'); % evaluate at these x values
   
end

f = f(:);

 
%
% find the nearest p points in N for each point in M
%
p = N - 1;  %integer between m+1 and N

d_matrix = zeros(N,1);

for i=1:M
    
    for j=1:N
        
        d_matrix(j) = 0;
        
        for k=1:m
            
            d_matrix(j) = d_matrix(j) + norm(y(i,k) - x(j,k));
               
        end
        
    end
    
    [z,index(i,:)] = sort(d_matrix);
    
    for j=1:p
        
        ip = (i-1)*p + j;
        
        points(ip,:) = x(index(i,j),:);
        
    end
    
end

%
% formulate the least square problem
%
for np = 1 : M
    
    A = [1 points((np-1)*p+1,:)];
    
    for i = (np-1)*p+2 : np*p
        
        A = [A; 1 points(i,:)];
        
    end
    
    B = f(index(np,1));
    
    for i=2:p
        
        B = [B; f(index(np,i))];
        
    end
    
    z = A \ B;
    
    if np == 1
        
       b_matrix = z(2:m+1);
       
    else
        
       b_matrix = [b_matrix z(2:m+1)];
       
    end
    
end

%construct the covariance matrix

C = 0;

for i=1 : M
    
    z = b_matrix(:,i);
    
    C = C + z * z';
    
end

C=C/M;

[W D] = eig(C);

[lambda_loclin, idx] = sort(diag(D), 'descend');

W = W(:,idx);


eta1 = W(:,1);

eta2 = W(:,2);



figure
semilogy(lambda_loclin./lambda_loclin(1), '-o');
set(gca, 'fontsize', 20);
title('local linear eigs');

%univariate SSP
figure
g1 = eta1'*x';
plot(g1, f, 'ko', 'markerfacecolor', 'k')
set(gca, 'fontsize', 20);
xlabel('y = <eta1, x>');
ylabel('f(x)');
title('local linear SSP');
