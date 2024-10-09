close all;
clear all;
clc; 


%% Physical parameters of the system:
M = 5;
m = 1;
g = 9.8;
l = 0.5;

t0 = 0;
T = 5;
dt = 0.001;
h = dt;
N = (T-t0)/dt; % descretização
k = 0:1:N;
t = t0 + k*dt;

%% Cost function array elements:
alpha_1 = 1;
alpha_2 = 1;
alpha_3 = 1e-6;

S = zeros(4,4);
P = zeros(4,4);
S(3,3)=alpha_1;
P(3,3)=alpha_2;

Q = alpha_3;

Y0 = [0;0;0.15;0]; 
y0 = Y0;


%% Problem linearization matrices, from the state ODE:
DF0 = [0,1,0,0; ...
       0, 0,-(g*m)/M,0; ...
       0,0,0,1; ...
       0,0,(g*M^2 + g*m*M)/(M^2*l),0];
F0 = [0;0;0;0];
G0 = [0;1/M;0;-1/(M*l)];

A = DF0;
B = G0;

M = [A, -B*Q^-1*B.'; -P, -A.'];

V0 = [eye(4); S];


%% RK4
% Butcher table:
aa = [0,0,0,0 ; 1/2,0,0,0 ; 0,1/2,0,0 ; 0,0,1,0];
bb = [1/6,1/3,1/3,1/6];
tt = [0 ; 1/2 ; 1/2 ; 1];
%%
V = zeros(8,4,N+1); % declaro o tensor do sistema que estou a resolver
for j = 1:1:4 % 4 colunas ou 4 sistemas de EDOs independentes
    v0 = V0(:,j);
    
    v = [];
    v = [v v0];
    
    for i = 1:N % varrimento no tempo

    vk = v(:,end);
    
    F1 = -M*vk;
    F2 = -M*( vk + aa(2,1)*h*F1);
    F3 = -M*( vk + aa(3,1)*h*F1+aa(3,2)*h*F2);
    F4 = -M*( vk + aa(4,1)*h*F1+aa(4,2)*h*F2+aa(4,3)*h*F3);
    
    vk_plus_1 = vk + h*(bb(1)*F1+bb(2)*F2+bb(3)*F3+bb(4)*F4);
    v = [v vk_plus_1];

    end 
    V(:,j,:) = flip(v,2);
    j
end


X = V(1:4,:,:);
Y = V(5:8,:,:);

R = zeros(4,4,N+1);
for i = 1:N+1
    R(:,:,i) = Y(:,:,i)*X(:,:,i)^-1;

end 

% maquina!!!!

y = [];
y = [y Y0];


for i = 1:N

yk = y(:,end);

F1 = (A-B*Q^-1*B.'*R(:,:,i))*( yk);
F2 = (A-B*Q^-1*B.'*R(:,:,i))*( yk + aa(2,1)*h*F1);
F3 = (A-B*Q^-1*B.'*R(:,:,i))*( yk + aa(3,1)*h*F1+aa(3,2)*h*F2);
F4 = (A-B*Q^-1*B.'*R(:,:,i))*( yk + aa(4,1)*h*F1+aa(4,2)*h*F2+aa(4,3)*h*F3);

yk_plus_1 = yk + h*(bb(1)*F1+bb(2)*F2+bb(3)*F3+bb(4)*F4);
y = [y yk_plus_1];
end
%% Control u(t):
u = zeros(N+1,1);
for i=1:N+1
    K = -Q^-1*B.'*R(:,:,i);
    u(i) = transpose(y(:,i))*transpose(K);
end 

%% Cost calculation J(u):
J = 0.5 * transpose(y0)*R(:,:,1)* y0
%J = 0.5 * transpose(y0)*R*y0
%% Plots
figure (1)
plot(t,y(1,:),t,y(2,:),t,y(3,:),t,y(4,:),'LineWidth',1.5)
legend('y1 - x','y2','y3 - \theta','y4')
grid;

figure (2)
plot(t,u,'LineWidth',1.5)
legend('u')
grid;



