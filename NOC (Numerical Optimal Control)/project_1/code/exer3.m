close all;
clear all;
clc; 

%% Physical parameters of the system
M = 5;
m = 1;
g = 9.8;
l = 0.5;


%% Discretization of time
t0 = 0;
T = 5;
h = 0.001;
dt =h;
N = (T-t0)/h; 
k = 0:1:N;
t = t0 + k*h;

step = 10^-18; % step of gradient descent

%% Cost function array elements:
alpha_1 = 1;
alpha_2 = 1;
alpha_3 = 1e-6;

S = zeros(4,4);
P = zeros(4,4);
S(3,3)=alpha_1;
P(3,3)=alpha_2;

Q = alpha_3;

%% Problem linearization matrices, from the state ODE:
% ODE Y'=F(Y,U)=AY+Bu arrays:
DF0 = [0,1,0,0; ...
       0, 0,-(g*m)/M,0; ...
       0,0,0,1; ...
       0,0,(g*M^2 + g*m*M)/(M^2*l),0];
F0 = [0;0;0;0];
G0 = [0;1/M;0;-1/(M*l)];

A = DF0;
B = G0;

%% Initializations
y = zeros(4,N+1); 
Y0 = [0;0;0.15;0];  % Initial State values
y(:,1) = Y0;

p = zeros(4,N+1);

 % Initial Control u(t):
u = ones(1,N+1)*0.3;


J_kplus1 = 10^9999;
%for k = 0:1:1
k=0;
while 1
    k
    
    % Calculate the state solution y(t), with Euler Implicit:
    y = zeros(4,N+1);
    y(:,1) = Y0;
    for n = 2:1:N+1
        y(:,n) = (eye(4)*(1/dt)-A)^-1*((1/dt)*y(:,n-1) + B*u(n));
    end 
    
    % Calculate the adjoint solution p(t), with Euler Implicit:
    p = zeros(4,N+1);
    p(:,end) = S*y(:,end);
    for n = N:-1:1
        p(:,n) = (-eye(4)*(1/dt)+A.')^-1*((-1/dt)*p(:,n+1) - P*y(:,n));
    end 

    % Gradient of Cost Function:
    grad_J = Q*u + B'*p;
    % Actualize the control u(t) with Gradient descent algorithm:
    u = u - step*grad_J;


    % Calculate the cost function J(u)
    f = diag(y'*P*y);
    q = diag(u'*Q*u);
 
    J_k = J_kplus1;
    J_kplus1 = 0.5*y(:,end)'*S*y(:,end)  + dt*0.5*(f(1)+f(end))+dt*sum(f(2:end-1)) + dt*0.5*(q(1)+q(end))+dt*sum(q(2:end-1))

    if abs(J_kplus1 - J_k) < 10^-7
        break;
    end 
    k = k+1;

end

%% Plots
figure (1)
plot(t,y(1,:),t,y(2,:),t,y(3,:),t,y(4,:),'LineWidth',1.5)
legend('y1 - x','y2','y3 - \theta','y4')
grid;

figure (2)
plot(t,u,'LineWidth',1.5)
legend('u')
grid;
