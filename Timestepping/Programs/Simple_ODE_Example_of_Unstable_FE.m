% A program to demonstrate instability of timestepping methods 
% when the timestep is inappropriately choosen.

%Differential equation: y'(t)=-y(t) y(t_0)=y_0
%Initial Condition, y(t_0)=1 where t_0=0
clear all; clc; clf;

%Grid
h=.1;
tmax=4;
Npoints = tmax/h;
lambda=.1;

%Initial Data
y0=1; 
t_0=0;
t(1)=t_0;
y_be(1)=y0;
y_fe(1)=y0;
y_imr(1)=y0;

for n=1:Npoints
	%Forward Euler
    y_fe(n+1)=y_fe(n)-lambda*h*y_fe(n);     
    	  %Backwards Euler
    y_be(n+1)=y_be(n)/(1+lambda*h);       
    	%Crank Nicolson
    y_imr(n+1)=y_imr(n)*(2-lambda*h)/(2+lambda*h) 
    t(n+1)=t(n)+h;
end

%Exact Solution
tt=[0:.001:tmax];
exact=exp(-lambda*tt);

%Plot
figure(1); clf; plot(tt,exact,'r-',t,y_fe,'b:',t,y_be,'g--',t,y_imr,'k-.');
xlabel time; ylabel y; 
legend('Exact','Forward Euler','Backward Euler',...
		'Crank Nicolson');
