function [ J,it_num] = LRR( data,alpha,lambda,mu,a,beta)
# min || X - XC ||_F^2 + lambda || C ||_*
#S via ADMM

if (nargin<3)
    mu = 1;
end

max_iterations =200;
n = size(data, 2);

J = zeros(n);
C = zeros(n);
Y = zeros(n);

tol1 = 2*1e-3;
tol2 = 1*1e-5;
T1 = C;
trIndex1 = double(T1 ~= 0);

var_X =data' * data;
var1 = norm((data - data * C),'fro');
var_R =inv(var1 + a * a);
A = inv(2 * var_X * var_R + mu * speye((n))+beta * speye((n)));

tol_1 = 1*10^-4;

for k = 1 : max_iterations
    
    Z_prev = J;
    J_prev = C;
    
    #W Y2
    [W,Y2] = Bound(alpha, beta, T1, trIndex1, tol1, tol2, max_iterations, 0, 1);
    
    #compute Ck+1
    C = A * (2 * var_X * var_R + mu * J + Y + beta * W + Y2);
    C = abs(C);
    #compute Jk+1
    V = C - (1/mu) * Y;
    [J, ~] = solve_nn(V, lambda / mu );   
    J(J<0) = 0;
    #compute Yk+1
    Y = Y + mu * (J - C);
    
    # Check convergence
    if (max(max(abs(C - J))) < tol_1&&max(max(abs(J - Z_prev))) < tol_1&&max(max(abs(C - J_prev))) < tol_1)
        break;
    end
    
end
it_num=k;
end