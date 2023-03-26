function X = RK4_InfiniteHorizon_States(X, X0, A, B, E, K)
X = zeros(6,1500);
X(:,1) = X0;
h = 0.01;
    for i = 1:1000-1
        f1 = LHS(X(:,i), A, B, E, K);
        f2 = LHS(X(:,i)+f1/2*h, A, B, E, K);
        f3 = LHS(X(:,i)+f2/2*h, A, B, E, K);
        f4 = LHS(X(:,i)+f3*h, A, B, E, K);
        X(:,i+1) = X(:,i) + h*(f1/6+(f2+f3)/3+f4/6);
    end
end

function dXdt = LHS(X, A, B, E, K)
    dXdt = E^(-1)*(A + B*K)*X;
end