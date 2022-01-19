function [T_recovery, T_recovery2] = Bound(alpha, beta, T, trIndex, tol1, tol2, maxiter, a, b)
X = T;
W = X;
Y = X;
i = 1;
stop1 = 1;
stop2 = 1;
while(stop1 > tol1 || stop2 > tol2)
    #the process of computing W
    tran = (1/beta) * (Y + alpha * (T.* trIndex)) + X;
    W = tran - (alpha/ (alpha + beta)) * (tran.* trIndex);
    W(W < a) = a; #a=0
    W(W > b) = b; #b=1 
    
    X_1 = svt(W- 1/beta* Y, 1/beta);
    
    Y = Y + beta * (X_1 - W);
    
    stop1_0 = stop1;
    stop1 = norm(X_1 - X, 'fro') / norm(X, 'fro');
    stop2 = abs(stop1 - stop1_0)/ max(1, abs(stop1_0));
    
    X = X_1;
    i = i+1;
   
    if i < maxiter
        iter = i - 1;
    else
        iter = maxiter;
        warning('reach maximum iteration~~do not converge!!!');
        break
    end
    
end

T_recovery = W;
T_recovery2 = Y;


end