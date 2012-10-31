% A program to solve the 2D nonlinear Schr\"{o}dinger equation using a
% splitting method 

clear all; format compact; format short;
set(0,'defaultaxesfontsize',30,'defaultaxeslinewidth',.7,...
    'defaultlinelinewidth',6,'defaultpatchlinewidth',3.7,'defaultaxesfontweight','bold')

% set up grid
tic
Lx = 20;        % period  2*pi*L
Ly = 20;        % period  2*pi*L
Nx = 2*256;     % number of harmonics
Ny = 2*256;     % number of harmonics
Nt = 100;       % number of time slices
dt = 5.0/Nt;    % time step

Es = 1.0;

% initialise variables
x = (2*pi/Nx)*(-Nx/2:Nx/2 -1)'*Lx;          % x coordinate
kx = 1i*[0:Nx/2-1 0 -Nx/2+1:-1]'/Lx;        % wave vector
y = (2*pi/Ny)*(-Ny/2:Ny/2 -1)'*Ly;          % y coordinate
ky = 1i*[0:Ny/2-1 0 -Ny/2+1:-1]'/Ly;        % wave vector
[xx,yy]=meshgrid(x,y);
[k2xm,k2ym]=meshgrid(kx.^2,ky.^2);
% initial conditions
u = exp(-(xx.^2+yy.^2));
v=fft2(u);
figure(1); clf; mesh(xx,yy,u); drawnow;
t=0; tdata(1)=t;

% mass
ma = fft2(abs(u).^2);
ma0 = ma(1,1);

% solve pde and plot results
for n =2:Nt+1
    vna=exp(0.5*1i*dt*(k2xm + k2ym)).*v;
    una=ifft2(vna);
    pot=Es*((abs(una)).^2);
    unb=exp(-1i*dt*pot).*una;
    vnb=fft2(unb);
    v=exp(0.5*1i*dt*(k2xm + k2ym)).*vnb;
    u=ifft2(v);
    t=(n-1)*dt;
    tdata(n)=t;
     if (mod(n,10)==0)
         figure(2); clf; mesh(xx,yy,abs(u).^2); title(num2str(t));
         drawnow;
         ma = fft2(abs(u).^2);
         ma = ma(1,1);
         test = log10(abs(1-ma/ma0))
     end
end
figure(4); clf; mesh(xx,yy,abs(u).^2);
toc