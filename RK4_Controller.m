function X = RK4_Controller(X, X0, A, B, E, R, Q)
X = zeros(size(A,1),size(A,2),size(A,3));
X(:,:, end) = X0;
h = 0.01;
    for i = size(X,3):-1:2
        f1 = RHS_Controller(X(:,:,i), A(:,:,i), B, E(:,:,i), R, Q);
        f2 = RHS_Controller(X(:,:,i)-f1/2*h, A(:,:,i), B, E(:,:,i), R, Q);
        f3 = RHS_Controller(X(:,:,i)-f2/2*h,  A(:,:,i), B, E(:,:,i), R, Q);
        f4 = RHS_Controller(X(:,:,i)-f3*h,  A(:,:,i), B, E(:,:,i), R, Q);
        X(:,:,i-1) = X(:,:,i) - h*(f1/6+(f2+f3)/3+f4/6);
    end
end

function dXdt = RHS_Controller(X, A, B, E, R, Q)
    dXdt = - E'^(-1)*A'*X - X*A*E^-1 + X*B*R^(-1)*B'*X - E'^(-1)*Q*E^-1;
end