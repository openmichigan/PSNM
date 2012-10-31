% A program to solve the 1D cubic Klein Gordon equation using a
% second order semi-explicit method 
% u_{tt}-u_{xx}+u=u^3
clear all; format compact; format short;
set(0,'defaultaxesfontsize',30,'defaultaxeslinewidth',.7,...
    'defaultlinelinewidth',6,'defaultpatchlinewidth',3.7,...
    'defaultaxesfontweight','bold')

% set up grid
tic
Lx = 64;            % period  2*pi*L
Nx = 4096;          % number of harmonics
Nt = 500;           % number of time slices
plotgap=10;         % time steps to take between plots
c=0.5;              % wave speed
dt = 5.00/Nt;       % time step

Es = 1.0;           % focusing (+1) or defocusing (-1) parameter
t=0; tdata(1)=t;

% initialise variables
x = (2*pi/Nx)*(-Nx/2:Nx/2 -1)'*Lx;          % x coordinate
kx = 1i*[0:Nx/2-1 0 -Nx/2+1:-1]'/Lx;        % wave vector

% initial conditions
u = sqrt(2)*sech((x-c*t)/sqrt(1-c^2));
uexact= sqrt(2)*sech((x-c*t)/sqrt(1-c^2));
uold=sqrt(2)*sech((x+c*dt)/sqrt(1-c^2));
v=fft(u,[],1);
vold=fft(uold,[],1);
figure(1); clf; 
% Plot data on
plot(x,u,'r+',x,uexact,'b-'); legend('numerical','exact');
title(num2str(t));  xlabel x; ylabel u;  drawnow;


% initial energy
vx=0.5*kx.*(v+vold);
ux=ifft(vx,[],1); 
Kineticenergy=0.5*abs( (u-uold)/dt).^2;           
Strainenergy=0.5*abs(ux).^2;           
Potentialenergy=0.5*abs(0.5*(u+uold)).^2 ...
                    -Es*0.25*((u+uold)*0.5).^4;           
Kineticenergy=fft(Kineticenergy,[],1);
Potentialenergy=fft(Potentialenergy,[],1);
Strainenergy=fft(Strainenergy,[],1);
EnKin(1)=Kineticenergy(1);
EnPot(1)=Potentialenergy(1);
EnStr(1)=Strainenergy(1);
En(1)=EnStr(1)+EnKin(1)+EnPot(1);
En0=En(1)

plotnum=1;
% solve pde and plot results

for n =1:Nt+1
    nonlin=u.^3;
    nonlinhat=fft(nonlin,[],1);
    vnew=(0.25*(kx.*kx -1).*(2*v+vold)...
                +(2*v-vold)/(dt*dt) +Es*nonlinhat)./...
            (1/(dt*dt) - (kx.*kx-1)*0.25 );
    unew=ifft(vnew,[],1);
    t=n*dt;
    if (mod(n,plotgap)==0)
        uexact=sqrt(2)*sech((x-c*t)/sqrt(1-c^2));
        figure(1); clf;
        plot(x,u,'r+',x,uexact,'b-'); legend('numerical','exact');
        title(num2str(t)); xlim([-6,6]); xlabel x; ylabel u;  drawnow;        
        tdata(plotnum+1)=t;
        vx=0.5*kx.*(v+vold);
        ux=ifft(vx,[],1); 
        Kineticenergy=0.5*abs( (u-uold)/dt).^2;           
        Strainenergy=0.5*abs(ux).^2;           
        Potentialenergy=0.5*abs(0.5*(u+uold)).^2 ...
                    -Es*0.25*((u+uold)*0.5).^4;           
        Kineticenergy=fft(Kineticenergy,[],1);
        Potentialenergy=fft(Potentialenergy,[],1);
        Strainenergy=fft(Strainenergy,[],1);
        EnKin(plotnum+1)=Kineticenergy(1);
        EnPot(plotnum+1)=Potentialenergy(1);
        EnStr(plotnum+1)=Strainenergy(1);
        En(plotnum+1)=EnStr(plotnum+1)+EnKin(plotnum+1)+EnPot(plotnum+1);
        Enchange(plotnum)=log(abs(1-En(1+plotnum)/En0));
        plotnum=plotnum+1;
    end
    % update old terms
    vold=v;
    v=vnew;
    uold=u;
    u=unew;
end
figure(4); clf; 
uexact=sqrt(2)*sech((x-c*t)/sqrt(1-c^2));
plot(x,u,'r+',x,uexact,'b-'); legend('numerical','exact');
title(num2str(t)); xlabel x; ylabel u;  drawnow;
max(abs(u-uexact))
figure(5); clf; plot(tdata,En,'r-',tdata,EnKin,'b:',tdata,EnPot,'g-.',tdata,EnStr,'y--'); 
xlabel time; ylabel Energy; legend('Total','Kinetic','Potential','Strain');
figure(6); clf; plot(tdata(2:end),Enchange,'r-'); xlabel time; ylabel('Energy change'); 

toc