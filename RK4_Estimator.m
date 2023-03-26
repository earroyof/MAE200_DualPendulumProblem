function P = RK4_Estimator(P, P0, A, C, E, Q1, Q2)
P = zeros(size(A,1),size(A,2),size(A,3));
P(:,:, 1) = P0;
h = 0.01;
    for i = 1:(size(A,3)-1)
        f1 = LHS_Estimator(P(:,:,i), A(:,:,i), C, E(:,:,i), Q1, Q2);
        f2 = LHS_Estimator(P(:,:,i)+f1/2*h, A(:,:,i), C, E(:,:,i), Q1, Q2);
        f3 = LHS_Estimator(P(:,:,i)+f2/2*h, A(:,:,i), C, E(:,:,i), Q1, Q2);
        f4 = LHS_Estimator(P(:,:,i)+f3*h, A(:,:,i), C, E(:,:,i), Q1, Q2);
        P(:,:,i+1) = P(:,:,i) + h*(f1/6+(f2+f3)/3+f4/6);
    end
end

function dPdt = LHS_Estimator(P, A, C, E, Q1, Q2)
    dPdt = E'^(-1)*A'*P + P*A*E^-1 - P*C'*Q2^(-1)*C*P + E'^(-1)*Q1*E^-1;
end