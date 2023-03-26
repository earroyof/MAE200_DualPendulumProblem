function X_hat = RK4_StateEstimatiomPerturbations(X_hat, X_prime, X, U, A, B, C, E, K, L)
X_hat = zeros(6,size(A,3));
X_hat(:, 1) = [1 1 1 1 1 1];
h = 0.01;
    for i = 1:(size(A,3)-1)
        f1 = LHS_StateEstimationPerturbations(X_hat(:,i), X_prime(:,i),X(1:6,i), U(i), A(:,:,i), B, C, E(:,:,i), K(:,:,i), L(:,:,i));
        f2 = LHS_StateEstimationPerturbations(X_hat(:,i)+f1/2*h, X_prime(:,i), X(1:6,i), U(i), A(:,:,i), B, C, E(:,:,i), K(:,:,i), L(:,:,i));
        f3 = LHS_StateEstimationPerturbations(X_hat(:,i)+f2/2*h, X_prime(:,i), X(1:6,i), U(i), A(:,:,i), B, C, E(:,:,i), K(:,:,i), L(:,:,i));
        f4 = LHS_StateEstimationPerturbations(X_hat(:,i)+f3*h, X_prime(:,i), X(1:6,i), U(i), A(:,:,i), B, C, E(:,:,i), K(:,:,i), L(:,:,i));
        X_hat(:,i+1) = X_hat(:,i) + h*(f1/6+(f2+f3)/3+f4/6);
    end
end

function dXhatDt = LHS_StateEstimationPerturbations(X_hat, X_prime, X, U, A, B, C, E, K, L)
    dXhatDt = (E^(-1)*A+L*C)*X_hat + (2*L*C-E^(-1)*B*K)*X + (L*C+E^(-1)*B*K)*X_prime + E^(-1)*B*U;
end