% A program to solve the u_t=u(u-1) using a
% splitting method

clear all; format compact; format short;
set(0,'defaultaxesfontsize',30,'defaultaxeslinewidth',.7,...
    'defaultlinelinewidth',6,'defaultpatchlinewidth',3.7,'defaultaxesfontweight','bold')
Nt = 1000;                          % number of time slices
tmax = 1;                           % maximum time
dt=tmax/Nt;                       % increment between times
time=(linspace(1,Nt,Nt)-1)*dt;      % time
uexact=4./(4+exp(time));% exact solution
u(1)=0.8

for i=1:Nt-1
    c=-1/u(i);
    utemp=-1/(c+dt);
    u(i+1)=utemp*exp(-dt);
end
figure(1)
plot(time,u,'r+',time,uexact,'b-');