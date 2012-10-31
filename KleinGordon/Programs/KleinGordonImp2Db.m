% A program to solve the 2D Klein Gordon equation using a
% second order implicit method 

clear all; format compact; format short;
set(0,'defaultaxesfontsize',30,'defaultaxeslinewidth',.7,...
    'defaultlinelinewidth',6,'defaultpatchlinewidth',3.7,...
    'defaultaxesfontweight','bold')

% set up grid
tic
Lx = 3;         % period  2*pi*L
Ly = 3;         % period  2*pi*L
Nx = 2*256;     % number of harmonics
Ny = 2*256;     % number of harmonics
Nt = 2000;        % number of time slices
dt = 50.0/Nt;  % time step
tol=0.1^(10);   % tolerance for fixed point iterations
plotgap=10;     % timesteps between plots

Es = 1.0; % focusing (+1) or defocusing (-1) parameter


% initialise variables
x = (2*pi/Nx)*(-Nx/2:Nx/2 -1)'*Lx;          % x coordinate
kx = 1i*[0:Nx/2-1 0 -Nx/2+1:-1]'/Lx;        % wave vector
y = (2*pi/Ny)*(-Ny/2:Ny/2 -1)'*Ly;          % y coordinate
ky = 1i*[0:Ny/2-1 0 -Ny/2+1:-1]'/Ly;        % wave vector
[xx,yy]=meshgrid(x,y);
[kxm,kym]=meshgrid(kx,ky);

% initial conditions
u = (0.5*exp(-(xx.^2+yy.^2))).*sin(10*xx+12*yy);
uold=u;
v=fft2(u); 
vold=fft2(uold); 
figure(1); clf; mesh(xx,yy,u); drawnow;
t=0; tdata(1)=t;

% initial energy
vx=0.5*kxm.*(v+vold);
vy=0.5*kym.*(v+vold);
ux=ifft2(vx); 
uy=ifft2(vy); 
ux=ifft2(vx); 
uy=ifft2(vy); 
Kineticenergy=0.5*abs( (u-uold)/dt).^2;           
Strainenergy=0.5*abs(ux).^2 +0.5*abs(uy).^2;           
Potentialenergy=0.5*abs(0.5*(u+uold)).^2 ...
                    -Es*0.25*((u+uold)*0.5).^4;           
Kineticenergy=fft2(Kineticenergy);
Potentialenergy=fft2(Potentialenergy);
Strainenergy=fft2(Strainenergy);
EnKin(1)=Kineticenergy(1,1);
EnPot(1)=Potentialenergy(1,1);
EnStr(1)=Strainenergy(1,1);
En(1)=EnStr(1)+EnKin(1)+EnPot(1);
En0=En(1)
plotnum=1;

% solve pde and plot results

for n =1:Nt+1
    nonlin=(u.^4 -uold.^4)./(u-uold+0.1^14);
    nonlinhat=fft2(nonlin);
    chg=1;
    unew=u;
    while (chg>tol)
        utemp=unew;
        vnew=(0.25*(kxm.^2 + kym.^2 -1).*(2*v+vold)...
             +(2*v-vold)/(dt*dt) +Es*nonlinhat)./...
            (1/(dt*dt) - (kxm.^2 +  kym.^2-1)*0.25 );
        unew=ifft2(vnew);
        nonlin=(unew.^4 -uold.^4)./(unew-uold+0.1^14);
        nonlinhat=fft2(nonlin);
        chg=max(abs(unew-utemp));
    end
    t=n*dt;
    if (mod(n,plotgap)==0)
        figure(1); clf; mesh(xx,yy,abs(u).^2);
        t
        tdata(plotnum+1)=t;
        vx=0.5*kxm.*(v+vold);
        vy=0.5*kym.*(v+vold);
        ux=ifft2(vx); 
        uy=ifft2(vy); 
        Kineticenergy=0.5*abs( (unew-u)/dt).^2;           
        Strainenergy=0.5*abs(ux).^2 +0.5*abs(uy).^2;           
        Potentialenergy=0.5*abs(0.5*(unew+u)).^2 ...
                    -Es*0.25*((unew+u)*0.5).^4;           
        Kineticenergy=fft2(Kineticenergy);
        Potentialenergy=fft2(Potentialenergy);
        Strainenergy=fft2(Strainenergy);
        EnKin(1+plotnum)=Kineticenergy(1,1);
        EnPot(1+plotnum)=Potentialenergy(1,1);
        EnStr(1+plotnum)=Strainenergy(1,1);
        En(1+plotnum)=EnStr(1+plotnum)+EnKin(1+plotnum)+EnPot(1+plotnum);
 
        Enchange(plotnum)=log(abs(1-En(1+plotnum)/En0));
        plotnum=plotnum+1;
    end
    % update old terms
    vold=v;
    v=vnew;
    uold=u;
    u=unew;
end
figure(5); clf; plot(tdata,En,'r-',tdata,EnKin,'b:',tdata,EnPot,'g-.',tdata,EnStr,'y--'); 
xlabel time; ylabel Energy; legend('Total','Kinetic','Potential','Strain');
figure(6); clf; plot(tdata(2:end),Enchange,'r-'); xlabel time; ylabel('Energy change'); 

figure(4); clf; mesh(xx,yy,abs(u).^2);
toc