close all;
clear all;
clc; 

%% Physical parameters of the system:
M = 5;
m = 1;
g = 9.8;
l = 0.5;

%% Cost function array elements:
alpha_1 = 1;
alpha_2 = 1.e-6;

P = eye(4);
P(1,1)=0;
P = P*alpha_1;

Q =alpha_2;

S = zeros(4,4);

%% Problem linearization matrices, from the state ODE:
DF0 = [0,1,0,0; ...
       0, 0,-(g*m)/M,0; ...
       0,0,0,1; ...
       0,0,(g*M^2 + g*m*M)/(M^2*l),0];
F0 = [0;0;0;0];
G0 = [0;1/M;0;-1/(M*l)];

A = DF0;
B = G0; 

%% Obtaining the matrix R, with the matlab sover
[K,R,L] = lqr(A,B,P,Q)
K = -Q^-1*B.'*R

%% Discretization of time
t0 = 0;
T = 10;
h = 0.0001;
N = (T-t0)/h; 
k = 0:1:N;
t = t0 + k*h;

%% RK4
% Butcher table:
aa = [0,0,0,0 ; 1/2,0,0,0 ; 0,1/2,0,0 ; 0,0,1,0];
bb = [1/6,1/3,1/3,1/6];
tt = [0 ; 1/2 ; 1/2 ; 1];

% Calculate the state solution y(t), with RK-4:
y0 = [0;0;0.15;0]; %Initial points
y = [];
y = [y y0];

for i = 1:N

yk = y(:,end);

F1 = F_y_ode_riccati(t(i)+tt(1)*h , yk,A,B,P,Q,R);
F2 = F_y_ode_riccati(t(i)+tt(2)*h , yk + aa(2,1)*h*F1,A,B,P,Q,R);
F3 = F_y_ode_riccati(t(i)+tt(3)*h , yk + aa(3,1)*h*F1+aa(3,2)*h*F2,A,B,P,Q,R);
F4 = F_y_ode_riccati(t(i)+tt(4)*h , yk + aa(4,1)*h*F1+aa(4,2)*h*F2+aa(4,3)*h*F3,A,B,P,Q,R);

yk_plus_1 = yk + h*(bb(1)*F1+bb(2)*F2+bb(3)*F3+bb(4)*F4);
y = [y yk_plus_1];
end

%% Calculation of control u(t):

%u= y(1,:)*K(1)+y(2,:)*K(2)+y(3,:)*K(3)+y(4,:)*K(4);
u = transpose(y)*transpose(K);

%% Cost calculation J(u):
J = 0.5 * transpose(y0)*R*y0

%% Plots
figure (1)
plot(t,y(1,:),t,y(2,:),t,y(3,:),t,y(4,:),'LineWidth',1.5)
legend('y1 - x','y2','y3 - \theta','y4')
grid;

figure (2)
plot(t,u,'LineWidth',1.5)
legend('u')
grid;


%% Auxiliar Functions
function F = F_y_ode_riccati(t,Y,A,B,P,Q,R)
   F = (A-B*Q^-1*B.'*R)*Y;
end
