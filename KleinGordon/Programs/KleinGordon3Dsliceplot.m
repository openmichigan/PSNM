% A program to solve the 3D Klein Gordon equation using a
% second order semi-explicit method 

clear all; format compact; format short;
set(0,'defaultaxesfontsize',30,'defaultaxeslinewidth',.7,...
    'defaultlinelinewidth',6,'defaultpatchlinewidth',3.7,...
    'defaultaxesfontweight','bold')

% set up grid
tic
Lx = 2;         % period  2*pi*L
Ly = 2;         % period  2*pi*L
Lz = 2;         % period  2*pi*L
Nx = 64;        % number of harmonics
Ny = 64;        % number of harmonics
Nz = 64;        % number of harmonics
Nt = 2000;       % number of time slices
plotgap=10;     
dt = 10.0/Nt;    % time step

Es = -1.0;  % focusing (+1) or defocusing (-1) parameter

% initialise variables
x = (2*pi/Nx)*(-Nx/2:Nx/2 -1)'*Lx;          % x coordinate
kx = 1i*[0:Nx/2-1 0 -Nx/2+1:-1]'/Lx;        % wave vector
y = (2*pi/Ny)*(-Ny/2:Ny/2 -1)'*Ly;          % y coordinate
ky = 1i*[0:Ny/2-1 0 -Ny/2+1:-1]'/Ly;        % wave vector
z = (2*pi/Nz)*(-Nz/2:Nz/2 -1)'*Lz;          % y coordinate
kz = 1i*[0:Nz/2-1 0 -Nz/2+1:-1]'/Lz;        % wave vector
[xx,yy,zz]=meshgrid(x,y,z);
[kxm,kym,kzm]=meshgrid(kx,ky,kz);

% initial conditions
u = 0.1*exp(-(xx.^2+(yy).^2+zz.^2));
uold=u;
v=fftn(u);
vold=v;
figure(1); clf; 
% coordinate slice to show plots on
sx=[0]; sy=[0]; sz=[-Lx*2*pi];
slice(xx,yy,zz,u,sx,sy,sz); colormap jet;
title(num2str(0)); colorbar('location','EastOutside'); drawnow;

xlabel('x'); ylabel('y'); zlabel('z'); 
axis equal; axis square;  view(3); drawnow;
t=0; tdata(1)=t;

% initial energy
vx=0.5*kxm.*(v+vold);
vy=0.5*kym.*(v+vold);      
vz=0.5*kzm.*(v+vold);      
ux=ifftn(vx); 
uy=ifftn(vy); 
uz=ifftn(vz); 
Kineticenergy=0.5*abs( (u-uold)/dt).^2;           
Strainenergy=0.5*abs(ux).^2 +0.5*abs(uy).^2+0.5*abs(uz).^2;           
Potentialenergy=0.5*abs(0.5*(u+uold)).^2 ...
                    -Es*0.25*((u+uold)*0.5).^4;           
Kineticenergy=fftn(Kineticenergy);
Potentialenergy=fftn(Potentialenergy);
Strainenergy=fftn(Strainenergy);
EnKin(1)=Kineticenergy(1,1);
EnPot(1)=Potentialenergy(1,1);
EnStr(1)=Strainenergy(1,1);
En(1)=EnStr(1)+EnKin(1)+EnPot(1);
En0=En(1)

plotnum=1;
% solve pde and plot results

for n =1:Nt+1
    nonlin=u.^3;
    nonlinhat=fftn(nonlin);
    vnew=(0.25*(kxm.^2 + kym.^2 + kzm.^2 -1).*(2*v+vold)...
                +(2*v-vold)/(dt*dt) +Es*nonlinhat)./...
            (1/(dt*dt) - (kxm.^2 + kzm.^2 + kym.^2 - 1)*0.25 );
    unew=ifftn(vnew);
    t=n*dt;
    if (mod(n,plotgap)==0)
        figure(1); clf; sx=[0]; sy=[0]; sz=[0];
        slice(xx,yy,zz,u,sx,sy,sz); colormap jet;
        title(num2str(t)); colorbar('location','EastOutside'); drawnow;
        xlabel('x'); ylabel('y'); zlabel('z'); 
        axis equal; axis square;  view(3); drawnow;
        tdata(plotnum+1)=t;
        t
        vx=0.5*kxm.*(v+vold);
        vy=0.5*kym.*(v+vold);      
        vz=0.5*kzm.*(v+vold);      
        ux=ifftn(vx); 
        uy=ifftn(vy); 
        uz=ifftn(vz); 
        Kineticenergy=0.5*abs( (u-uold)/dt).^2;           
        Strainenergy=0.5*abs(ux).^2 +0.5*abs(uy).^2+0.5*abs(uz).^2;           
        Potentialenergy=0.5*abs(0.5*(u+uold)).^2 ...
                    -Es*0.25*((u+uold)*0.5).^4;           
        Kineticenergy=fftn(Kineticenergy);
        Potentialenergy=fftn(Potentialenergy);
        Strainenergy=fftn(Strainenergy);
        EnKin(plotnum+1)=Kineticenergy(1,1,1);
        EnPot(plotnum+1)=Potentialenergy(1,1,1);
        EnStr(plotnum+1)=Strainenergy(1,1,1);
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
% coordinate slice to show plots on
sx=[0]; sy=[0]; sz=[0];
slice(xx,yy,zz,u,sx,sy,sz); colormap jet;
title(num2str(t)); colorbar('location','EastOutside'); drawnow;

xlabel('x'); ylabel('y'); zlabel('z'); 
axis equal; axis square;  view(3); drawnow;

figure(5); clf; plot(tdata,En,'r-',tdata,EnKin,'b:',tdata,EnPot,'g-.',tdata,EnStr,'y--'); 
xlabel time; ylabel Energy; legend('Total','Kinetic','Potential','Strain');
figure(6); clf; plot(tdata(2:end),Enchange,'r-'); xlabel time; ylabel('Energy change'); 

toc