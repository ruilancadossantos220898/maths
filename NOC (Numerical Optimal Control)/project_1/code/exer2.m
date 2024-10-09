close all;
close;
clear all;
clear;
clc; 

%% Physical parameters of the system:
M = 5;
m = 1;
g = 9.8;
l = 0.5;

Y0 = [0;0;0.15;0]; % Initial State
Y = [];
Y = [Y0 Y]

t0 = 0;
tN = 10;
h = 0.001;
N = (tN-t0)/h; % discretization
k = 0:1:N;
t = t0 + k*h;


%% Runge-Kutta 4th order
% Butcher table:
aa = [0,0,0,0 ; 1/2,0,0,0 ; 0,1/2,0,0 ; 0,0,1,0];
bb = [1/6,1/3,1/3,1/6];
tt = [0 ; 1/2 ; 1/2 ; 1];

for i = 1:N
    Yk = Y(:,end);
    
    F1 = F_exer2(t(i)+tt(1)*h , Yk);
    F2 = F_exer2(t(i)+tt(2)*h , Yk + aa(2,1)*h*F1);
    F3 = F_exer2(t(i)+tt(3)*h , Yk + aa(3,1)*h*F1+aa(3,2)*h*F2);
    F4 = F_exer2(t(i)+tt(4)*h , Yk + aa(4,1)*h*F1+aa(4,2)*h*F2+aa(4,3)*h*F3);
    
    Yk_plus_1 = Yk + h*(bb(1)*F1+bb(2)*F2+bb(3)*F3+bb(4)*F4);
        
    Y = [Y Yk_plus_1];
end


%% Plots 

figure (1)
plot(t,Y(1,:),t,Y(2,:),t,Y(3,:),t,Y(4,:),'LineWidth',1.5)
legend('y1 - x','y2','y3 - \theta','y4')
grid


%% Auxiliar Functions

function F = F_exer2(t,Y)
   y1 = Y(1);
   y2 = Y(2);
   y3 = Y(3);
   y4 = Y(4);
   u = 0 ;
   M = 5;
   m = 1;
   g = 9.8;
   l = 0.5;

   DF0 = [0,1,0,0; ...
          0, 0,-(g*m)/M,0; ...
          0,0,0,1; ...
          0,0,(g*M^2 + g*m*M)/(M^2*l),0];

   F0 = [0;0;0;0];
   G0 = [0;1/M;0;1/(M*l)];

   F = DF0*Y + F0 + G0*u;
end

