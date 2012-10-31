% A program to solve the 3D nonlinear Schr\"{o}dinger equation using a
% splitting method 

clear all; format compact; format short;
set(0,'defaultaxesfontsize',30,'defaultaxeslinewidth',.7,...
    'defaultlinelinewidth',6,'defaultpatchlinewidth',3.7,...
    'defaultaxesfontweight','bold')

% set up grid
tic
Lx = 4;         % period  2*pi*L
Ly = 4;         % period  2*pi*L
Lz = 4;         % period  2*pi*L
Nx = 64;        % number of harmonics
Ny = 64;        % number of harmonics
Nz = 64;        % number of harmonics
Nt = 100;       % number of time slices
dt = 1.0/Nt;    % time step

Es = 1.0;  % focusing or defocusing parameter

% initialise variables
x = (2*pi/Nx)*(-Nx/2:Nx/2 -1)'*Lx;          % x coordinate
kx = 1i*[0:Nx/2-1 0 -Nx/2+1:-1]'/Lx;        % wave vector
y = (2*pi/Ny)*(-Ny/2:Ny/2 -1)'*Ly;          % y coordinate
ky = 1i*[0:Ny/2-1 0 -Ny/2+1:-1]'/Ly;        % wave vector
z = (2*pi/Nz)*(-Nz/2:Nz/2 -1)'*Lz;          % y coordinate
kz = 1i*[0:Nz/2-1 0 -Nz/2+1:-1]'/Lz;        % wave vector
[xx,yy,zz]=meshgrid(x,y,z);
[k2xm,k2ym,k2zm]=meshgrid(kx.^2,ky.^2,kz.^2);

% initial conditions
u = exp(-(xx.^2+yy.^2+zz.^2));
v=fftn(u);
figure(1); clf; UP = abs(u).^2;
p1 = patch(isosurface(x,y,z,UP,.0025),...
    'FaceColor','yellow','EdgeColor','none');
p2 = patch(isocaps(x,y,z,UP,.0025),...
    'FaceColor','interp','EdgeColor','none');
isonormals(UP,p1); lighting phong;
xlabel('x'); ylabel('y'); zlabel('z'); 
axis equal; axis square; view(3); drawnow;
t=0; tdata(1)=t;

% mass
ma = fftn(abs(u).^2);
ma0 = ma(1,1,1);

% solve pde and plot results

for n =2:Nt+1
    vna=exp(0.5*1i*dt*(k2xm + k2ym + k2zm)).*v;
    una=ifftn(vna);
    pot=Es*((abs(una)).^2);
    unb=exp(-1i*dt*pot).*una;
    vnb=fftn(unb);
    v=exp(0.5*1i*dt*(k2xm + k2ym + k2zm)).*vnb;
    u=ifftn(v);
    t=(n-1)*dt;
    tdata(n)=t;
    if (mod(n,10)==0)
        figure(1); clf; UP = abs(u).^2;
        p1 = patch(isosurface(x,y,z,UP,.0025),...
            'FaceColor','yellow','EdgeColor','none');
        p2 = patch(isocaps(x,y,z,UP,.0025),...
            'FaceColor','interp','EdgeColor','none');
        isonormals(UP,p1); lighting phong;
        xlabel('x'); ylabel('y'); zlabel('z'); 
        axis equal; axis square; view(3); drawnow;
        ma = fftn(abs(u).^2);
        ma = ma(1,1,1);  test = log10(abs(1-ma/ma0))
     end
end
figure(4); clf; UP = abs(u).^2;
p1 = patch(isosurface(x,y,z,UP,.0025),...
    'FaceColor','yellow','EdgeColor','none');
p2 = patch(isocaps(x,y,z,UP,.0025),...
    'FaceColor','interp','EdgeColor','none');
isonormals(UP,p1); lighting phong;
xlabel('x'); ylabel('y'); zlabel('z'); 
axis equal; axis square;  view(3); drawnow;
toc