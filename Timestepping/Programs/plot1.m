% A program to demonstrate instability of timestepping methods 
% when the timestep is inappropriately choosen.

%Differential equation: y'(t)=-y(t) y(t_0)=y_0
%Initial Condition, y(t_0)=1 where t_0=0
clear all; clc; clf;
set(0,'defaultaxesfontsize',25,'defaultaxeslinewidth',.7,...
    'defaultlinelinewidth',6,'defaultpatchlinewidth',3.7,...
    'defaultaxesfontweight','bold')
%Grid
h=1.0;
tmax=4;
Npoints = tmax/h;
lambda=-1;

%Initial Data
y0=1; 
t_0=0;
t1(1)=t_0;
y_fe1(1)=y0;

for n=1:Npoints
	%Forward Euler
    y_fe1(n+1)=y_fe1(n)-lambda*h*y_fe1(n);     
    t1(n+1)=t1(n)+h;
end

%Grid
h=0.10;
tmax=4;
Npoints = tmax/h;
lambda=-1;

%Initial Data
y0=1; 
t_0=0;
t2(1)=t_0;
y_fe2(1)=y0;

for n=1:Npoints
	%Forward Euler
    y_fe2(n+1)=y_fe2(n)-lambda*h*y_fe2(n);     
    t2(n+1)=t2(n)+h;
end

%Exact Solution
tt=[0:.001:tmax];
exact=exp(-lambda*tt);

%Plot
figure(1); clf; plot(tt,exact,'r-',t1,y_fe1,'b:',t2,y_fe2,'g--');
xlabel Time; ylabel y(t); 
legend('Exact','h=1','h=0.1');
