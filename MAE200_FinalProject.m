%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Phase A 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Constants
T = 3; %sec

s.h=0.01; s.N=T/s.h; s.mc=10; s.t=[0:s.N]*s.h;               % STEP 0: initialize simulation
s.m1=1; s.L1=1;    s.ell1=s.L1; s.I1=s.m1*s.ell1^2/3; % system, & derived parameters
s.m2=0.5; s.L2=0.5;  s.ell2=s.L2; s.I2=s.m2*s.ell2^2/3; alpha=0.1; 
s.B=[0; 0; 0; 1; 0; 0]; s.Q=diag([1 1 1 1 1 1]); s.R = alpha^2; s.QT=diag([5 40 10 .1 60 10]); %s.R=0; %[5 40 25 .1 55 15]

[u_k,x_k] = Dual_Inverted_Pendulum(T,s);

% state(:,1) = x_k(1:6,end);
% for i = 1:5
%     [u_k, x_k] = dual_inverted_pendulum(T,u_k);
%     state(:,i+1) = x_k(1:6,end);
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Compute A and E
A = zeros(6,6,size(u_k,1));
E = zeros(6,6,size(u_k,1));
for i = 1:length(x_k)
    A(:,:,i) = Compute_A(x_k(:,i),s);
    E(:,:,i) = Compute_E(x_k(:,i), s);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% March the Controller Riccati Eq backwards to get X & Solve for controller
% gain K
X = zeros(6,6,size(u_k,1));
K = zeros(1,6,size(u_k,1));
X0=eye(6);
% EQ 22.13b
X = RK4_Controller(X, X0, A, s.B, E, s.R, s.Q);
for i = 1:length(x_k) 
    K(:,:,i) = -s.R^(-1)*s.B'*X(:,:,i);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% March the Estimator Riccati Eq backwards to get P & Solve for observer
% gain L
P = zeros(6,6,size(u_k,1));
L = zeros(6,6,size(u_k,1));
P0 = eye(6);
Q1 = eye(6);
Q2 = eye(6);
C = diag([1 1 1 0 0 0]);
% Eq 22.30
P = RK4_Estimator(P, P0, A, C, E, Q1, Q2);
for i = 1:length(x_k)
    L(:,:,i) = -P(:,:,i)*C'*Q2^-1;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plotting x_k + hx' & x_k to check how good the K is 
X_prime = zeros(6,size(u_k,1));
X_prime = RK4_StatePerturbations(X_prime, x_k, u_k, A, s.B, E, K);

plotAxis = 1:1:size(u_k,1);
% Simulated states = nominal states + perturbations

figure 
subplot(2,3,1)
plot(plotAxis,x_k(1,:),plotAxis,x_k(1,:) + s.h*X_prime(1,:))
title('State 1')
legend('Open Loop State','Control Closed Loop State','Location','SouthWest')

subplot(2,3,2)
plot(plotAxis,x_k(2,:),plotAxis,x_k(2,:) + s.h*X_prime(2,:))
title('State 2')
legend('Open Loop State','Control Closed Loop State','Location','SouthWest')

subplot(2,3,3)
plot(plotAxis,x_k(3,:),plotAxis,x_k(3,:) + s.h*X_prime(3,:))
title('State 3')
legend('Open Loop State','Control Closed Loop State','Location','SouthWest')

subplot(2,3,4)
plot(plotAxis,x_k(4,:),plotAxis,x_k(4,:) + s.h*X_prime(4,:))
title('State 4')
legend('Open Loop State','Control Closed Loop State','Location','SouthWest')

subplot(2,3,5)
plot(plotAxis,x_k(5,:),plotAxis,x_k(5,:) + s.h*X_prime(5,:))
title('State 5')
legend('Open Loop State','Control Closed Loop State','Location','SouthWest')

subplot(2,3,6)
plot(plotAxis,x_k(6,:),plotAxis,x_k(6,:) + s.h*X_prime(6,:))
title('State 6')
legend('Open Loop State','Control Closed Loop State','Location','SouthWest')

sgtitle('Phase A: LQR')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plotting x_k + hx_hat & x_k to check how good the L is 
X_hat = zeros(6,size(u_k,1));
X_hat = RK4_StateEstimatiomPerturbations(X_hat, X_prime, x_k, u_k, A, s.B, C, E, K, L);

figure 
subplot(2,3,1)
plot(plotAxis,x_k(1,:),plotAxis,x_k(1,:) + s.h*X_hat(1,:))
title('State 1: x')
legend('Open Loop State','Control + Estimation State','Location','SouthWest')

subplot(2,3,2)
plot(plotAxis,x_k(2,:),plotAxis,x_k(2,:) + s.h*X_hat(2,:))
title('State 2: theta1')
legend('Open Loop State','Control + Estimation State','Location','SouthWest')

subplot(2,3,3)
plot(plotAxis,x_k(3,:),plotAxis,x_k(3,:) + s.h*X_hat(3,:))
title('State 3: theta2')
legend('Open Loop State','Control + Estimation State','Location','SouthWest')

subplot(2,3,4)
plot(plotAxis,x_k(4,:),plotAxis,x_k(4,:) + s.h*X_hat(4,:))
title('State 4: dx')
legend('Open Loop State','Control + Estimation State','Location','SouthWest')

subplot(2,3,5)
plot(plotAxis,x_k(5,:),plotAxis,x_k(5,:) + s.h*X_hat(5,:))
title('State 5: dtheta1')
legend('Open Loop State','Control + Estimation State','Location','SouthWest')

subplot(2,3,6)
plot(plotAxis,x_k(6,:),plotAxis,x_k(6,:) + s.h*X_hat(6,:))
title('State 6: dtheta2')
legend('Open Loop State','Control + Estimation State','Location','SouthWest')

sgtitle('Phase A: LQG')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Phase B - Infinite Horizon Problem
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
g=9.8;
I = eye(3); 
Z = zeros(3);

% New E and A matrices for linearized model
E_inf = [I Z; Z [s.mc+s.m1+s.m2         -s.m1*s.ell1 -s.m2*s.ell2;
           -s.m1*s.ell1  s.I1+s.m1*s.ell1^2             0        ;
           -s.m2*s.ell2          0           s.I2+s.m2*s.ell2^2]];
A_inf=[Z I; [0 0 0 ; 0 s.m1*g*s.ell1 0; 0  0 s.m2*g*s.ell2] Z];

[X_icare, Kx_icare, Lx_icare] = icare(E_inf^(-1)*A_inf, E_inf^(-1)*s.B, s.Q, s.R);
K_inf = -s.R^(-1)*(E_inf^(-1)*s.B)'*X_icare;

% Get the states
X_inf = zeros(6,1500);
X_inf = RK4_InfiniteHorizon_States(X_inf, x_k(1:6,end), A_inf, s.B, E_inf, K_inf);

plotAxis_inf = 1:1:size(X_inf,2);
figure 
subplot(2,3,1)
plot(plotAxis_inf,X_inf(1,:))
title('State 1: x')

subplot(2,3,2)
plot(plotAxis_inf,X_inf(2,:))
title('State 2: theta1')

subplot(2,3,3)
plot(plotAxis_inf,X_inf(3,:))
title('State 3: theta2')

subplot(2,3,4)
plot(plotAxis_inf,X_inf(4,:))
title('State 4: dx')

subplot(2,3,5)
plot(plotAxis_inf,X_inf(5,:))
title('State 5: dtheta1')

subplot(2,3,6)
plot(plotAxis_inf,X_inf(6,:))
title('State 6: dtheta2')

sgtitle('Phase B: LQR')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[P_icare, Kp_icare, Lp_icare] = icare((E_inf^(-1)*A_inf)', C, Q1, Q2);
L_inf = -P_icare*C'*Q2^-1;

% Get the state estimations
X_hat_inf = zeros(6,1500);
X_hat_inf = RK4_InfiniteHorizon_StateEstimations(X_hat_inf, X_inf, x_k(1:6,end), A_inf, s.B, C, K_inf, L_inf);

figure 
subplot(2,3,1)
plot(plotAxis_inf,X_hat_inf(1,:))
title('State 1: x')

subplot(2,3,2)
plot(plotAxis_inf,X_hat_inf(2,:))
title('State 2: theta1')

subplot(2,3,3)
plot(plotAxis_inf,X_hat_inf(3,:))
title('State 3: theta2')

subplot(2,3,4)
plot(plotAxis_inf,X_hat_inf(4,:))
title('State 4: dx')

subplot(2,3,5)
plot(plotAxis_inf,X_hat_inf(5,:))
title('State 5: dtheta1')

subplot(2,3,6)
plot(plotAxis_inf,X_hat_inf(6,:))
title('State 6: dtheta2')

sgtitle('Phase B: LQG')