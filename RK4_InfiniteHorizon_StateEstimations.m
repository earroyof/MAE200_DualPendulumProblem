function X_hat_inf = RK4_InfiniteHorizon_StateEstimations(X_hat_inf, X, X0, A, B, C, K, L)
X_hat_inf = zeros(6,1500);
X_hat_inf(:,1) = X0;
h = 0.01;
    for i = 1:1000-1
        f1 = LHS(X_hat_inf(:,i), X(:,i), A, B, C, K, L);
        f2 = LHS(X_hat_inf(:,i)+f1/2*h, X(:,i), A, B, C, K, L);
        f3 = LHS(X_hat_inf(:,i)+f2/2*h, X(:,i), A, B, C, K, L);
        f4 = LHS(X_hat_inf(:,i)+f3*h, X(:,i), A, B, C, K, L);
        X_hat_inf(:,i+1) = X_hat_inf(:,i) + h*(f1/6+(f2+f3)/3+f4/6);
    end
end

function dXhatDt = LHS(X_hat, X, A, B, C, K, L)
    dXhatDt = (A + L*C)*X_hat + (B*K - L*C)*X;
end