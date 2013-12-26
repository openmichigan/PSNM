% A program to solve the 3D incompressible magnetohydrodynamic equations 
% with periodic boundary
% conditions. The program is based on the Orszag-Patterson algorithm as
% documented on pg. 98 of C. Canuto, M.Y. Hussaini, A. Quarteroni and
% T.A. Zhang "Spectral Methods: Evolution to Complex Geometries and 
% Applications to Fluid Dynamics" Springer (2007)
%
% The Helmholtz decomposition is used to project the magnetic field onto a
% divergence free subspace. Initial work on this has been done with Damian
% San Roman Alerigi
%
% For plotting, the program uses vol3d, which is available at:
%
% http://www.mathworks.com/matlabcentral/fileexchange/22940-vol3d-v2
%
clear all; format compact; format short;
set(0,'defaultaxesfontsize',30,'defaultaxeslinewidth',.7,...
    'defaultlinelinewidth',6,'defaultpatchlinewidth',3.7,...
    'defaultaxesfontweight','bold')

% set up grid
tic
Lx = 1;         % period  2*pi*L
Ly = 1;         % period  2*pi*L
Lz = 1;         % period  2*pi*L
Nx = 32;        % number of harmonics
Ny = 32;        % number of harmonics
Nz = 32;        % number of harmonics
Nt = 50;        % number of time slices
dt = 2/Nt;    % time step
t=0;            % initial time
Re = 100.0;       % Reynolds number
Rem = 100.0;      % Magnetic Reynolds number
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

% initial condition
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

bx=sin(sl*yy).*sin(sm.*zz);
by=sin(sm.*zz);
bz=cos(sk*xx).*cos(sl*yy);

bxhat=fftn(bx);
byhat=fftn(by);
bzhat=fftn(bz);

bxx=ifftn(bxhat.*kxm);bxy=ifftn(bxhat.*kym);bxz=ifftn(bxhat.*kzm);
byx=ifftn(byhat.*kxm);byy=ifftn(byhat.*kym);byz=ifftn(byhat.*kzm);
bzx=ifftn(bzhat.*kxm);bzy=ifftn(bzhat.*kym);bzz=ifftn(bzhat.*kzm);

% calculate fluid vorticity for plotting
omegax=wy-vz; omegay=uz-wx; omegaz=vx-uy; 
omegatot=omegax.^2+omegay.^2+omegaz.^2;
% calculate magnetic vorticity for plotting    
omegabx=bxy-byz; omegaby=bxz-bzx; omegabz=byx-bxy; 
omegabtot=omegabx.^2+omegaby.^2+omegabz.^2;
figure(1); clf; n=0;
subplot(2,1,1); title(['|Fluid omega|^2 ',num2str(t)]);
H = vol3d('CData',omegatot,'texture','3D','XData',x,'YData',y,'ZData',z);
xlabel('x'); ylabel('y'); zlabel('z'); colorbar;
axis equal; axis square; view(3); 
    
subplot(2,1,2); title(['|Magnetic omega|^2 ',num2str(t)]);
H = vol3d('CData',omegabtot,'texture','3D','XData',x,'YData',y,'ZData',z);
xlabel('x'); ylabel('y'); zlabel('z'); colorbar;
axis equal; axis square; view(3); 


for n=1:Nt
    uold=u; uxold=ux; uyold=uy; uzold=uz; 
    vold=v; vxold=vx; vyold=vy; vzold=vz; 
    wold=w; wxold=wx; wyold=wy; wzold=wz; 

    bxold=bx; bxxold=bxx; bxyold=bxy; bxzold=bxz;
    byold=by; byxold=byx; byyold=byy; byzold=byz;
    bzold=bz; bzxold=bzx; bzyold=bzy; bzzold=bzz;
    
    rhsuhatfix=(1/dt+(0.5/Re)*(k2xm+k2ym+k2zm)).*uhat;
    rhsvhatfix=(1/dt+(0.5/Re)*(k2xm+k2ym+k2zm)).*vhat;
    rhswhatfix=(1/dt+(0.5/Re)*(k2xm+k2ym+k2zm)).*what;
    
    rhsbxhatfix=(1/dt+(0.5/Rem)*(k2xm+k2ym+k2zm)).*bxhat;
    rhsbyhatfix=(1/dt+(0.5/Rem)*(k2xm+k2ym+k2zm)).*byhat;
    rhsbzhatfix=(1/dt+(0.5/Rem)*(k2xm+k2ym+k2zm)).*bzhat;

    chg=1; t=t+dt;
    while (chg>tol)
        nonlinu=0.25*((u+uold).*(ux+uxold)...
                      +(v+vold).*(uy+uyold)...
                      +(w+wold).*(uz+uzold)...
                      -(bx+bxold).*(bxx+bxxold)...
                      -(by+byold).*(bxy+bxyold)...
                      -(bz+bzold).*(bxz+bxzold));
        nonlinv=0.25*((u+uold).*(vx+vxold)...
                      +(v+vold).*(vy+vyold)...
                      +(w+wold).*(vz+vzold)...
                      -(bx+bxold).*(byx+byxold)...
                      -(by+byold).*(byy+byyold)...
                      -(bz+bzold).*(byz+byzold));
        nonlinw=0.25*((u+uold).*(wx+wxold)...
                      +(v+vold).*(wy+wyold)...
                      +(w+wold).*(wz+wzold)...
                      -(bx+bxold).*(bzx+bzxold)...
                      -(by+byold).*(bzy+bzyold)...
                      -(bz+bzold).*(bzz+bzzold));
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
        
        nonlinu=0.25*((u+uold).*(bxx+bxxold)...
                      +(v+vold).*(bxy+bxyold)...
                      +(w+wold).*(bxz+bxzold)...
                      -(bx+bxold).*(ux+uxold)...
                      -(by+byold).*(uy+uyold)...
                      -(bz+bzold).*(uz+uzold));
        nonlinv=0.25*((u+uold).*(byx+byxold)...
                      +(v+vold).*(byy+byyold)...
                      +(w+wold).*(byz+byzold)...
                      -(bx+bxold).*(vx+vxold)...
                      -(by+byold).*(vy+vyold)...
                      -(bz+bzold).*(vz+vzold));
        nonlinw=0.25*((u+uold).*(bzx+bzxold)...
                      +(v+vold).*(bzy+bzyold)...
                      +(w+wold).*(bzz+bzzold)...
                      -(bx+bxold).*(wx+wxold)...
                      -(by+byold).*(wy+wyold)...
                      -(bz+bzold).*(wz+wzold));
        nonlinuhat=fftn(nonlinu);
        nonlinvhat=fftn(nonlinv);
        nonlinwhat=fftn(nonlinw);
        phat=-1.0*(kxm.*nonlinuhat+kym.*nonlinvhat+kzm.*nonlinwhat)...
            ./(k2xm+k2ym+k2zm+0.1^13);
        bxhat=(rhsbxhatfix-nonlinuhat-kxm.*phat)...
            ./(1/dt - (0.5/Re)*(k2xm+k2ym+k2zm));
        byhat=(rhsbyhatfix-nonlinvhat-kym.*phat)...
            ./(1/dt - (0.5/Re)*(k2xm+k2ym+k2zm));
        bzhat=(rhsbzhatfix-nonlinwhat-kzm.*phat)...
            ./(1/dt - (0.5/Re)*(k2xm+k2ym+k2zm));
        bxx=ifftn(bxhat.*kxm);bxy=ifftn(bxhat.*kym);bxz=ifftn(bxhat.*kzm);
        byx=ifftn(byhat.*kxm);byy=ifftn(byhat.*kym);byz=ifftn(byhat.*kzm);
        bzx=ifftn(bzhat.*kxm);bzy=ifftn(bzhat.*kym);bzz=ifftn(bzhat.*kzm);
        bxtemp=bx; bytemp=by; bztemp=bz;
        bx=ifftn(bxhat); by=ifftn(byhat); bz=ifftn(bzhat);
        
        chg=max(abs(utemp-u))   +max(abs(vtemp-v))  +max(abs(wtemp-w))+...
            max(abs(bxtemp-bx)) +max(abs(bytemp-by))+max(abs(bztemp-bz));
    end
    % calculate vorticity for plotting
    omegax=wy-vz; omegay=uz-wx; omegaz=vx-uy; 
    omegatot=omegax.^2+omegay.^2+omegaz.^2;
    
    % calculate magnetic vorticity for plotting
    omegabx=bxy-byz; omegaby=bxz-bzx; omegabz=byx-bxy; 
    omegabtot=omegabx.^2+omegaby.^2+omegabz.^2;
    figure(1); clf; 
    subplot(2,1,1); title(['|Fluid omega|^2 ',num2str(t)]);
    H = vol3d('CData',omegatot,'texture','3D','XData',x,'YData',y,'ZData',z);
    xlabel('x'); ylabel('y'); zlabel('z'); colorbar;
    axis equal; axis square; view(3); 
      
    subplot(2,1,2); title(['|Magnetic omega|^2 ',num2str(t)]);
    H = vol3d('CData',omegabtot,'texture','3D','XData',x,'YData',y,'ZData',z);
    xlabel('x'); ylabel('y'); zlabel('z'); colorbar;
    axis equal; axis square; view(3); 
       
end



