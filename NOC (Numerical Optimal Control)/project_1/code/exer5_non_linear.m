close all;
clear all;
clc; 
%% RK4 for ODE's
% Butcher table:
aa = [0,0,0,0 ; 1/2,0,0,0 ; 0,1/2,0,0 ; 0,0,1,0];
bb = [1/6,1/3,1/3,1/6];
tt = [0 ; 1/2 ; 1/2 ; 1];
%% Physic parameters of the problem:
M = 5;
m = 1;
g = 9.8;
l = 0.5;

t0 = 0;
T = 1;
step = 10^-2;

dt = 0.001;
h = dt;
N = (T-t0)/dt; % discretization
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

y = zeros(4,N+1);
Y0 = [0;0;0.15;0]; 
y(:,1) = Y0;

p = zeros(4,N+1);



% Initial Control:
u = ones(1,N+1)*0;

k=0;
J_kplus1 = 10^200;
%for k = 0:1:100000
while 1
    

    % Calculate the state solution y(t), with RK-4:
    y = zeros(4,N+1);
    y(:,1) = Y0;
    for i = 2:N+1 % RK-4 for ODE IVP:
    yk = y(:,i-1);
    
    F1 = function_y_ode(t(i)+tt(1)*h , yk,u(i));
    F2 = function_y_ode(t(i)+tt(2)*h , yk + aa(2,1)*h*F1,u(i));
    F3 = function_y_ode(t(i)+tt(3)*h , yk + aa(3,1)*h*F1+aa(3,2)*h*F2,u(i));
    F4 = function_y_ode(t(i)+tt(4)*h , yk + aa(4,1)*h*F1+aa(4,2)*h*F2+aa(4,3)*h*F3,u(i));
    
    yk_plus_1 = yk + h*(bb(1)*F1+bb(2)*F2+bb(3)*F3+bb(4)*F4);
    y(:,i) = yk_plus_1;
    end    

    % Calculate the adjoint solution p(t), with RK-4:
    pz = zeros(4,N+1);
    pz(:,1) = S*y(:,end);
    for i = 2:N+1
    
       pzk = pz(:,i-1);
       
       F1 = -function_p_ode(T-t(i)-tt(1)*h , pzk, y(:,N+2-i) , u(N+2-i), P);
       F2 = -function_p_ode(T-t(i)-tt(2)*h , pzk + aa(2,1)*h*F1, y(:,N+2-i) , u(N+2-i), P);
       F3 = -function_p_ode(T-t(i)-tt(3)*h , pzk + aa(3,1)*h*F1+aa(3,2)*h*F2,y(:,N+2-i) , u(N+2-i), P);
       F4 = -function_p_ode(T-t(i)-tt(4)*h , pzk + aa(4,1)*h*F1+aa(4,2)*h*F2+aa(4,3)*h*F3, y(:,N+2-i) , u(N+2-i), P);
       
       pzk_plus_1 = pzk + h*(bb(1)*F1+bb(2)*F2+bb(3)*F3+bb(4)*F4);
       pz(:,i) = pzk_plus_1;
    end
    p = flip(pz,2);

    % Gradient of Cost Function:
    grad_J = zeros(1,N+1);
    for i = 1:1:N+1
        grad_J(i) = Q*u(i) + grad_u_of_f(y(:,i))'*p(:,i);
    end 
    % Actualize the control u(t) with Gradient descent algorithm:
    u = u - step*grad_J;
    %norm_u = norm(u)
    %norm_grad_J = norm(grad_J)
    %k
    % Calculate the cost function J(u)
    J_k = J_kplus1;
    J_kplus1 = calculate_cost_function(y, u, S, P, Q, dt)
    if abs(J_kplus1 - J_k) < 10^-7
        break;
    end 
    if J_kplus1 > J_k 
        %step = step/2
        %fprintf('step too big')
        %pause 
    end 
figure(1)
subplot(2,1,1);
plot(t,u,'LineWidth',1.5)
legend('u')
grid
%ylim([0 100])
subplot(2,1,2); 
plot(t,y(1,:),t,y(2,:),t,y(3,:),t,y(4,:),'LineWidth',1.5)
legend('y1 - x','y2','y3 - \theta','y4')
grid
%ylim([-1 5])
     


k = k+1;
end
k
%% Plots 

figure (2)
plot(t,y(1,:),t,y(2,:),t,y(3,:),t,y(4,:),'LineWidth',1.5)
legend('y1 - x','y2','y3 - \theta','y4')
grid;
figure (3)
plot(t,u,'LineWidth',1.5)
legend('u')
grid;


%% Auxiliar Functions
function F = function_y_ode(t,Y,u)
   M = 5;
   m = 1;
   g = 9.8;
   l = 0.5;
   y1 = Y(1);
   y2 = Y(2);
   y3 = Y(3);
   y4 = Y(4);
F = [y2; ...
     (l*m*sin(y3)*abs(y4)^2 + u)/(- m*cos(y3)^2 + M + m) - (g*m*cos(y3)*sin(y3))/(- m*cos(y3)^2 + M + m); ...
     y4; ...
     (g*sin(y3)*(M + m))/(- l*m*cos(y3)^2 + M*l + l*m) - (cos(y3)*(l*m*sin(y3)*abs(y4)^2 + u))/(- l*m*cos(y3)^2 + M*l + l*m)];
end 



function F = function_p_ode(t,p,Y,u,P)
   M = 5;
   m = 1;
   g = 9.8;
   l = 0.5;
   y1 = Y(1);
   y2 = Y(2);
   y3 = Y(3);
   y4 = Y(4);

 
   grad_y_of_f = [0, 1,                                                                                                                                                                                                                                                                                                                                               0,                                                                       0; ...
                  0, 0,                                                                     (g*m*sin(y3)^2)/(- m*cos(y3)^2 + M + m) - (g*m*cos(y3)^2)/(- m*cos(y3)^2 + M + m) - (2*m*cos(y3)*sin(y3)*(l*m*sin(y3)*abs(y4)^2 + u))/(- m*cos(y3)^2 + M + m)^2 + (2*g*m^2*cos(y3)^2*sin(y3)^2)/(- m*cos(y3)^2 + M + m)^2 + (l*m*abs(y4)^2*cos(y3))/(- m*cos(y3)^2 + M + m),                (2*l*m*abs(y4)*sign(y4)*sin(y3))/(- m*cos(y3)^2 + M + m); ...
                  0, 0,                                                                                                                                                                                                                                                                                                                                               0,                                                                       1; ...
                  0, 0, (sin(y3)*(l*m*sin(y3)*abs(y4)^2 + u))/(- l*m*cos(y3)^2 + M*l + l*m) + (g*cos(y3)*(M + m))/(- l*m*cos(y3)^2 + M*l + l*m) - (l*m*abs(y4)^2*cos(y3)^2)/(- l*m*cos(y3)^2 + M*l + l*m) + (2*l*m*cos(y3)^2*sin(y3)*(l*m*sin(y3)*abs(y4)^2 + u))/(- l*m*cos(y3)^2 + M*l + l*m)^2 - (2*g*l*m*cos(y3)*sin(y3)^2*(M + m))/(- l*m*cos(y3)^2 + M*l + l*m)^2, -(2*l*m*abs(y4)*cos(y3)*sign(y4)*sin(y3))/(- l*m*cos(y3)^2 + M*l + l*m)];
        
   grad_y_of_l = P*Y;

   F = -grad_y_of_f.'*p - grad_y_of_l;
end 

function F = grad_u_of_f(Y)
   M = 5;
   m = 1;
   g = 9.8;
   l = 0.5;
   y1 = Y(1);
   y2 = Y(2);
   y3 = Y(3);
   y4 = Y(4);

F =[
 
                                                      0;
                              1/(- m*cos(y3)^2 + M + m);
                                                      0;
                 -cos(y3)/(- l*m*cos(y3)^2 + M*l + l*m);];
end




function J = calculate_cost_function(y, u, S, P, Q, dt)
    f = diag(y'*P*y);
    q = diag(u'*Q*u);
    J = 0.5*y(:,end)'*S*y(:,end)  + dt*0.5*(f(1)+f(end))+dt*sum(f(2:end-1)) + dt*0.5*(q(1)+q(end))+dt*sum(q(2:end-1));

end 