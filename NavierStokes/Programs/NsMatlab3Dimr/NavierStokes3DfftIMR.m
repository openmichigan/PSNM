% A program to solve the 3D Navier stokes equations with periodic boundary
% conditions. The program is based on the Orszag-Patterson algorithm as
% documented on pg. 98 of C. Canuto, M.Y. Hussaini, A. Quarteroni and
% T.A. Zhang "Spectral Methods: Evolution to Complex Geometries and 
% Applications to Fluid Dynamics" Springer (2007)
%
% The exact solution used to check the numerical method is in 
% A. Shapiro "The use of an exact solution of the Navier-Stokes equations 
% in a validation test of a three-dimensional nonhydrostatic numerical 
% model" Monthly Weather Review vol. 121 pp. 2420-2425 (1993) 

clear all; format compact; format short;
set(0,'defaultaxesfontsize',30,'defaultaxeslinewidth',.7,...
    'defaultlinelinewidth',6,'defaultpatchlinewidth',3.7,...
    'defaultaxesfontweight','bold')

% set up grid
tic
Lx = 1;         % period  2*pi*L
Ly = 1;         % period  2*pi*L
Lz = 1;         % period  2*pi*L
Nx = 64;        % number of harmonics
Ny = 64;        % number of harmonics
Nz = 64;        % number of harmonics
Nt = 10;       % number of time slices
dt = 0.2/Nt;    % time step
t=0;            % initial time
Re = 1.0;  % Reynolds number
tol=10^(-10);
% initialise variables
x = (2*pi/Nx)*(-Nx/2:Nx/2 -1)'*Lx;          % x coordinate
kx = 1i*[0:Nx/2-1 0 -Nx/2+1:-1]'/Lx;        % wave vector
y = (2*pi/Ny)*(-Ny/2:Ny/2 -1)'*Ly;          % y coordinate
ky = 1i*[0:Ny/2-1 0 -Ny/2+1:-1]'/Ly;        % wave vector
z = (2*pi/Nz)*(-Nz/2:Nz/2 -1)'*Lz;          % y coordinate
kz = 1i*[0:Nz/2-1 0 -Nz/2+1:-1]'/Lz;        % wave vector
[xx,yy,zz]=meshgrid(x,y,z);
[kxm,kym,kzm]=meshgrid(kx,ky,kz);
[k2xm,k2ym,k2zm]=meshgrid(kx.^2,ky.^2,kz.^2);

% initial conditions for Taylor-Green vortex
% theta=0;
% u=(2/sqrt(3))*sin(theta+2*pi/3)*sin(xx).*cos(yy).*cos(zz);
% v=(2/sqrt(3))*sin(theta-2*pi/3)*cos(xx).*sin(yy).*cos(zz);
% w=(2/sqrt(3))*sin(theta)*cos(xx).*cos(yy).*sin(zz);

% exact solution
sl=1; sk=1; sm=1; lamlkm=sqrt(sl.^2+sk.^2+sm.^2);
u=-0.5*(lamlkm*sl*cos(sk*xx).*sin(sl*yy).*sin(sm.*zz)...
            +sm*sk*sin(sk*xx).*cos(sl*yy).*cos(sm.*zz))...
            .*exp(-t*(lamlkm^2)/Re);
        
v=0.5*(lamlkm*sk*sin(sk*xx).*cos(sl*yy).*sin(sm.*zz)...
            -sm*sl*cos(sk*xx).*sin(sl*yy).*cos(sm.*zz))...
            *exp(-t*(lamlkm^2)/Re);
      
w=cos(sk*xx).*cos(sl*yy).*sin(sm*zz)*exp(-t*(lamlkm^2)/Re);

uhat=fftn(u);
vhat=fftn(v);
what=fftn(w);

ux=ifftn(uhat.*kxm);uy=ifftn(uhat.*kym);uz=ifftn(uhat.*kzm);
vx=ifftn(vhat.*kxm);vy=ifftn(vhat.*kym);vz=ifftn(vhat.*kzm);
wx=ifftn(what.*kxm);wy=ifftn(what.*kym);wz=ifftn(what.*kzm);

% calculate vorticity for plotting
omegax=wy-vz; omegay=uz-wx; omegaz=vx-uy; 
omegatot=omegax.^2+omegay.^2+omegaz.^2;
figure(1); clf; n=0;
subplot(2,2,1); title(['omega x ',num2str(n*dt)]);
p1 = patch(isosurface(x,y,z,omegax,.0025),...
            'FaceColor','interp','EdgeColor','none','FaceAlpha',0.3);
p2 = patch(isocaps(x,y,z,omegax,.0025),...
            'FaceColor','interp','EdgeColor','none','FaceAlpha',0.1);
        isonormals(omegax,p1); lighting phong;
xlabel('x'); ylabel('y'); zlabel('z'); 
axis equal; axis square; view(3); colorbar;
subplot(2,2,2); title(['omega y ',num2str(n*dt)]);
p1 = patch(isosurface(x,y,z,omegay,.0025),...
            'FaceColor','interp','EdgeColor','none','FaceAlpha',0.3);
p2 = patch(isocaps(x,y,z,omegay,.0025),...
            'FaceColor','interp','EdgeColor','none','FaceAlpha',0.1);
        isonormals(omegay,p1); lighting phong;
xlabel('x'); ylabel('y'); zlabel('z'); 
axis equal; axis square; view(3); colorbar;
subplot(2,2,3); title(['omega z ',num2str(n*dt)]);
p1 = patch(isosurface(x,y,z,omegaz,.0025),...
            'FaceColor','interp','EdgeColor','none','FaceAlpha',0.3);
p2 = patch(isocaps(x,y,z,omegaz,.0025),...
            'FaceColor','interp','EdgeColor','none','FaceAlpha',0.1);
        isonormals(omegaz,p1); lighting phong;
xlabel('x'); ylabel('y'); zlabel('z'); 
axis equal; axis square; view(3); colorbar;
subplot(2,2,4); title(['|omega|^2 ',num2str(n*dt)]);
p1 = patch(isosurface(x,y,z,omegatot,.0025),...
            'FaceColor','interp','EdgeColor','none','FaceAlpha',0.3);
p2 = patch(isocaps(x,y,z,omegatot,.0025),...
            'FaceColor','interp','EdgeColor','none','FaceAlpha',0.1);
        isonormals(omegatot,p1); lighting phong;
xlabel('x'); ylabel('y'); zlabel('z'); colorbar;
axis equal; axis square; view(3);


for n=1:Nt
    uold=u; uxold=ux; uyold=uy; uzold=uz; 
    vold=v; vxold=vx; vyold=vy; vzold=vz; 
    wold=w; wxold=wx; wyold=wy; wzold=wz; 
    rhsuhatfix=(1/dt+(0.5/Re)*(k2xm+k2ym+k2zm)).*uhat;
    rhsvhatfix=(1/dt+(0.5/Re)*(k2xm+k2ym+k2zm)).*vhat;
    rhswhatfix=(1/dt+(0.5/Re)*(k2xm+k2ym+k2zm)).*what;
    chg=1; t=t+dt;
    while (chg>tol)
        nonlinu=0.25*((u+uold).*(ux+uxold)...
                      +(v+vold).*(uy+uyold)...
                      +(w+wold).*(uz+uzold));
        nonlinv=0.25*((u+uold).*(vx+vxold)...
                      +(v+vold).*(vy+vyold)...
                      +(w+wold).*(vz+vzold));
        nonlinw=0.25*((u+uold).*(wx+wxold)...
                      +(v+vold).*(wy+wyold)...
                      +(w+wold).*(wz+wzold));
        nonlinuhat=fftn(nonlinu);
        nonlinvhat=fftn(nonlinv);
        nonlinwhat=fftn(nonlinw);
        phat=-1.0*(kxm.*nonlinuhat+kym.*nonlinvhat+kzm.*nonlinwhat)...
            ./(k2xm+k2ym+k2zm+0.1^13);
        uhat=(rhsuhatfix-nonlinuhat-kxm.*phat)...
            ./(1/dt - (0.5/Re)*(k2xm+k2ym+k2zm));
        vhat=(rhsvhatfix-nonlinvhat-kym.*phat)...
            ./(1/dt - (0.5/Re)*(k2xm+k2ym+k2zm));
        what=(rhswhatfix-nonlinwhat-kzm.*phat)...
            ./(1/dt - (0.5/Re)*(k2xm+k2ym+k2zm));
        ux=ifftn(uhat.*kxm);uy=ifftn(uhat.*kym);uz=ifftn(uhat.*kzm);
        vx=ifftn(vhat.*kxm);vy=ifftn(vhat.*kym);vz=ifftn(vhat.*kzm);
        wx=ifftn(what.*kxm);wy=ifftn(what.*kym);wz=ifftn(what.*kzm);
        utemp=u; vtemp=v; wtemp=w;
        u=ifftn(uhat); v=ifftn(vhat); w=ifftn(what);
        chg=max(abs(utemp-u))+max(abs(vtemp-v))+max(abs(wtemp-w));
    end
    % calculate vorticity for plotting
    omegax=wy-vz; omegay=uz-wx; omegaz=vx-uy; 
    omegatot=omegax.^2+omegay.^2+omegaz.^2;
    figure(1); clf; 
    subplot(2,2,1); title(['omega x ',num2str(t)]);
    p1 = patch(isosurface(x,y,z,omegax,.0025),...
            'FaceColor','interp','EdgeColor','none','FaceAlpha',0.3);
    p2 = patch(isocaps(x,y,z,omegax,.0025),...
            'FaceColor','interp','EdgeColor','none','FaceAlpha',0.1);
        isonormals(omegax,p1); lighting phong;
    xlabel('x'); ylabel('y'); zlabel('z'); 
    axis equal; axis square; view(3); colorbar;
    subplot(2,2,2); title(['omega y ',num2str(t)]);
    p1 = patch(isosurface(x,y,z,omegay,.0025),...
            'FaceColor','interp','EdgeColor','none','FaceAlpha',0.3);
    p2 = patch(isocaps(x,y,z,omegay,.0025),...
            'FaceColor','interp','EdgeColor','none','FaceAlpha',0.1);
        isonormals(omegay,p1); lighting phong;
    xlabel('x'); ylabel('y'); zlabel('z'); 
    axis equal; axis square; view(3); colorbar;
    subplot(2,2,3); title(['omega z ',num2str(t)]);
    p1 = patch(isosurface(x,y,z,omegaz,.0025),...
            'FaceColor','interp','EdgeColor','none','FaceAlpha',0.3);
    p2 = patch(isocaps(x,y,z,omegaz,.0025),...
            'FaceColor','interp','EdgeColor','none','FaceAlpha',0.1);
        isonormals(omegaz,p1); lighting phong;
    xlabel('x'); ylabel('y'); zlabel('z'); 
    axis equal; axis square; view(3); colorbar;
    subplot(2,2,4); title(['|omega|^2 ',num2str(t)]);
    p1 = patch(isosurface(x,y,z,omegatot,.0025),...
            'FaceColor','interp','EdgeColor','none','FaceAlpha',0.3);
    p2 = patch(isocaps(x,y,z,omegatot,.0025),...
            'FaceColor','interp','EdgeColor','none','FaceAlpha',0.1);
        isonormals(omegatot,p1); lighting phong;
    xlabel('x'); ylabel('y'); zlabel('z'); colorbar;
    axis equal; axis square; view(3);
end
toc

uexact=-0.5*(lamlkm*sl*cos(sk*xx).*sin(sl*yy).*sin(sm.*zz)...
            +sm*sk*sin(sk*xx).*cos(sl*yy).*cos(sm.*zz))...
            .*exp(-t*(lamlkm^2)/Re);
        
vexact=0.5*(lamlkm*sk*sin(sk*xx).*cos(sl*yy).*sin(sm.*zz)...
            -sm*sl*cos(sk*xx).*sin(sl*yy).*cos(sm.*zz))...
            *exp(-t*(lamlkm^2)/Re);
      
wexact=cos(sk*xx).*cos(sl*yy).*sin(sm*zz)*exp(-t*(lamlkm^2)/Re);


error=  max(max(max(abs(u-uexact))))+...
        max(max(max(abs(v-vexact))))+...
        max(max(max(abs(w-wexact))))
