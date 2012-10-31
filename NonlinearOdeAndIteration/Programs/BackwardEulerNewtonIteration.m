% A program to solve y'=y^2 using the backward Euler
% method and Newton iteration
% This is not optimized and is very simple

set(0,'defaultaxesfontsize',25,'defaultaxeslinewidth',.7,...
'defaultlinelinewidth',6,'defaultpatchlinewidth',3.7,...
    'defaultaxesfontweight','bold')

n=100000; % Number of timesteps
Tmax=0.99; % Maximum time
y0=1; % Initial value
tol=0.1^10; % Tolerance for fixed point iterations
dt=Tmax/n; % Time step
y=zeros(1,n); % vector for discrete solution
t=zeros(1,n); % vector for times of discrete solution
y(1)=y0;
t(1)=0;
tic,        % start timing
for i=1:n
    yold=y(i); ynew=y(i); err=1;
    while err>tol
        ynew=yold-(yold-y(i)-dt*yold^2)/(1-2*dt*yold);
        err=abs(ynew-yold);
        yold=ynew;
    end
    y(i+1)=ynew;
    t(i+1)=t(i)+dt;
end
toc,        % stop timing
yexact=1./(1-t); max(abs(y-yexact)), % print maximum error
figure(1); clf; plot(t,y,'r+',t,yexact,'b-.');
xlabel Time; ylabel Solution; 
legend('Backward Euler','Exact solution');
title('Numerical solution of dy/dt=y^2');