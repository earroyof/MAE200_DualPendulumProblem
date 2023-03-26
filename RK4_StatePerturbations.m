function X_prime = RK4_StatePerturbations(X_prime, X, U, A, B, E, K)
X_prime = zeros(6,size(A,3));
X_prime(:, 1) = [1 1 1 1 1 1];
h = 0.01;
    for i = 1:(size(A,3)-1)
        f1 = LHS_StatePerturbations(X_prime(:,i), X(1:6,i), U(i), A(:,:,i), B, E(:,:,i), K(:,:,i));
        f2 = LHS_StatePerturbations(X_prime(:,i)+f1/2*h, X(1:6,i), U(i), A(:,:,i), B, E(:,:,i), K(:,:,i));
        f3 = LHS_StatePerturbations(X_prime(:,i)+f2/2*h, X(1:6,i), U(i), A(:,:,i), B, E(:,:,i), K(:,:,i));
        f4 = LHS_StatePerturbations(X_prime(:,i)+f3*h, X(1:6,i), U(i), A(:,:,i), B, E(:,:,i), K(:,:,i));
        X_prime(:,i+1) = X_prime(:,i) + h*(f1/6+(f2+f3)/3+f4/6);
    end
end

function dXprimeDt = LHS_StatePerturbations(X_prime, X, U, A, B, E, K)
    dXprimeDt = -E^(-1)*B*K*X + E^(-1)*B*U + E^(-1)*(A+B*K)*X_prime;
end