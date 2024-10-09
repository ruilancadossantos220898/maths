close all;
clear all;
clc;

%% Physical parameters of the system
M = 5;
m = 1;
g = 9.8;
l = 0.5;

Y0 = [0;0;0.15;0];
Y = [];
Y = [Y0 Y]

t0 = 0;
T = 10;
h = 0.01;
N = (T-t0)/h; % discretization
k = 0:1:N;
t = t0 + k*h;

%% RK4
% Butcher table:
aa = [0,0,0,0 ; 1/2,0,0,0 ; 0,1/2,0,0 ; 0,0,1,0];
bb = [1/6,1/3,1/3,1/6];
tt = [0 ; 1/2 ; 1/2 ; 1];

for i = 1:N
    Yk = Y(:,end);
    
    F1 = F_exer1(t(i)+tt(1)*h , Yk);
    F2 = F_exer1(t(i)+tt(2)*h , Yk + aa(2,1)*h*F1);
    F3 = F_exer1(t(i)+tt(3)*h , Yk + aa(3,1)*h*F1+aa(3,2)*h*F2);
    F4 = F_exer1(t(i)+tt(4)*h , Yk + aa(4,1)*h*F1+aa(4,2)*h*F2+aa(4,3)*h*F3);
    
    Yk_plus_1 = Yk + h*(bb(1)*F1+bb(2)*F2+bb(3)*F3+bb(4)*F4);
    Y = [Y Yk_plus_1];
end

%% Plots

figure (1)
plot(t,Y(1,:),t,Y(2,:),t,Y(3,:),t,Y(4,:),'LineWidth',1.5)
grid
legend('y1 - x','y2','y3 - \theta','y4')


%% Auxiliar Functions

function F = F_exer1(t,Y)
   y1 = Y(1);
   y2 = Y(2);
   y3 = Y(3);
   y4 = Y(4);

   M = 5;
   m = 1;
   g = 9.8;
   l = 0.5;
   
   f1 = y2;
   f2 = (l*m*abs(y4)^2*sin(y3))/(- m*cos(y3)^2 + M + m) - (g*m*cos(y3)*sin(y3))/(- m*cos(y3)^2 + M + m);
   f3 = y4;
   f4 = (g*sin(y3)*(M + m))/(- l*m*cos(y3)^2 + M*l + l*m) - (l*m*abs(y4)^2*cos(y3)*sin(y3))/(- l*m*cos(y3)^2 + M*l + l*m) ;

   F = [f1;f2;f3;f4];
end

