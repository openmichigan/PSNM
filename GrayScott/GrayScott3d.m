% A program to solve the Gray-Scott equations using
% splitting method. The nonlinear equations are solved using the implicit
% midpoint rule
clear all; format compact, format short,
set(0,'defaultaxesfontsize',17,'defaultaxeslinewidth',.7,...
    'defaultlinelinewidth',3.5,'defaultpatchlinewidth',3.5)
 
Nx =64; % number of modes
Ny=64;
Nz=64;
% set up grid
Nt=4;
tmax=0.04;
Lx = 1.5;       % period  2*pi * L
Ly =1.5;
Lz=1.5;
dt = tmax/Nt;   % number of time slices
tol=0.1^12;
A = 0.04;
B = 0.1;
Du = 1.;
Dv = 1.;
% initialise variables
x = (2*pi/Nx)*(-Nx/2:Nx/2 -1)'*Lx;      % x coordinate
y = (2*pi/Ny)*(-Nx/2:Ny/2 -1)'*Ly;
z = (2*pi/Nz)*(-Nx/2:Nz/2 -1)'*Lz;
kx = i*[0:Nx/2-1 0 -Nx/2+1:-1]'/Lx;     % wave vector
ky = i*[0:Ny/2-1 0 -Ny/2+1:-1]'/Ly;
kz = i*[0:Nz/2-1 0 -Nz/2+1:-1]'/Lz;
[xx,yy,zz]=meshgrid(x,y,z);
[kxx,kyy,kzz]=meshgrid(kx,ky,kz);
% initial conditions
t=0; tdata(1)=t;
u = 0.2 + exp(-2*(xx.^2+yy.^2+zz.^2));
v = 0.1 + exp(-4*(xx.^2+yy.^2+zz.^2-0.01).^2);
 
gamma=[1];
Ahat=A*fftn(ones(size(xx)));
t=0;
figure(11); clf;
subplot(2,1,1);
H = vol3d('CData',real(u),'texture','3D','XData',x,'YData',y,'ZData',z);
xlabel('x'); ylabel('y'); zlabel('z'); colorbar;
axis equal; axis square; view(3);
 
xlabel('x'); ylabel('y'); zlabel('z');
axis equal; axis square; view(3); drawnow;
title(['Time ',num2str(t)]);
subplot(2,1,2);
H = vol3d('CData',real(v),'texture','3D','XData',x,'YData',y,'ZData',z);
xlabel('x'); ylabel('y'); zlabel('z'); colorbar;
axis equal; axis square; view(3);
 
xlabel('x'); ylabel('y'); zlabel('z');
axis equal; axis square; view(3); drawnow;
filename=['./pictures1/',num2str(10000000),'.jpg'];
saveas(11,filename)
 
 
for n=1:Nt
    chg=1;
    for m=1:1
        % use fixed point iteration to solve nonlinear system in real space
        chg=1;
        uold=u; vold=v;
        while (chg>tol)
            utemp=u; vtemp=v;
            umean=0.5*(u+uold);
            vmean=0.5*(v+vold);
            u=uold+0.5*gamma(m)*dt*(-umean.*vmean.^2);
            v=vold+0.5*gamma(m)*dt*umean.*vmean.^2;
            chg=max(abs(u-utemp))+max(abs(v-vtemp));
        end
        uhat=fftn(u); vhat=fftn(v); % solve linear part exactly in Fourier space
        uhat=exp(gamma(m)*dt*(-A+Du*(kxx.^2+kyy.^2+kzz.^2))).*...
            (uhat-Ahat./(A+Du*(kxx.^2+kyy.^2+kzz)))+Ahat./...
            (A+Du*(kxx.^2+kyy.^2+kzz.^2));
        vhat=exp(gamma(m)*dt*(-B+Dv*(kxx.^2+kyy.^2+kzz.^2))).*vhat;
        u=ifftn(uhat); v=ifftn(vhat);
        % use fixed point iteration to solve nonlinear system in real space
        chg=1;
        uold=u; vold=v;
        while (chg>tol)
            utemp=u; vtemp=v;
            umean=0.5*(u+uold);
            vmean=0.5*(v+vold);
            u=uold+0.5*gamma(m)*dt*(-umean.*vmean.^2);
            v=vold+0.5*gamma(m)*dt*umean.*vmean.^2;
            chg=max(abs(u-utemp))+max(abs(v-vtemp));
        end
    end
    t=n*dt;
    figure(11); clf;
    subplot(2,1,1);
    H = vol3d('CData',real(u),'texture','3D','XData',x,'YData',y,'ZData',z);
    xlabel('x'); ylabel('y'); zlabel('z'); colorbar;
    axis equal; axis square; view(3);
 
    xlabel('x'); ylabel('y'); zlabel('z');
    axis equal; axis square; view(3); drawnow;
    title(['Time ',num2str(t)]);
    subplot(2,1,2);
    H = vol3d('CData',real(v),'texture','3D','XData',x,'YData',y,'ZData',z);
    xlabel('x'); ylabel('y'); zlabel('z'); colorbar;
    axis equal; axis square; view(3);
 
    xlabel('x'); ylabel('y'); zlabel('z');
    axis equal; axis square; view(3); drawnow;
    filename=['./pictures1/',num2str(1000000+n),'.jpg'];
    saveas(11,filename)
 
end
