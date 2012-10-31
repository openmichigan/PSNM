% Numerical solution of the 2D incompressible Navier-Stokes on a
% Square Domain [0,1]x[0,1] using a Fourier pseudo-spectral method
% and Crank-Nicolson timestepping. The numerical solution is compared to
% the exact Taylor-Green Vortex solution of the Navier-Stokes equations
%
%Periodic free-slip boundary conditions and Initial conditions:
    %u(x,y,0)=sin(2*pi*x)cos(2*pi*y)
    %v(x,y,0)=-cos(2*pi*x)sin(2*pi*y)
%Analytical Solution:
    %u(x,y,t)=sin(2*pi*x)cos(2*pi*y)exp(-8*pi^2*t/Re)
    %v(x,y,t)=-cos(2*pi*x)sin(2*pi*y)exp(-8*pi^2*t/Re)
clear all; format compact; format short; clc; clf;

Re=1;%Reynolds number

%grid
Nx=64; h=1/Nx; x=h*(1:Nx);
Ny=64; h=1/Ny; y=h*(1:Ny)';
[xx,yy]=meshgrid(x,y);

%initial conditions
u=sin(2*pi*xx).*cos(2*pi*yy);
v=-cos(2*pi*xx).*sin(2*pi*yy);
u_y=-2*pi*sin(2*pi*xx).*sin(2*pi*yy);
v_x=2*pi*sin(2*pi*xx).*sin(2*pi*yy);
omega=v_x-u_y;

dt=0.0025; t(1)=0; tmax=.1;
nplots=ceil(tmax/dt);

%wave numbers for derivatives
k_x=2*pi*(1i*[(0:Nx/2-1)  0 1-Nx/2:-1]');
k_y=2*pi*(1i*[(0:Ny/2-1)  0 1-Ny/2:-1]);
k2x=k_x.^2;
k2y=k_y.^2;

%wave number grid for multiplying matricies
[kxx,kyy]=meshgrid(k2x,k2y);
[kx,ky]=meshgrid(k_x,k_y);

% use a high tolerance so time stepping errors
% are not dominated by errors in solution to nonlinear
% system
tol=10^(-10);

%compute \hat{\omega}^{n+1,k}
omegahat=fft2(omega);
%nonlinear term
nonlinhat=fft2(u.*ifft2(omegahat.*kx)+v.*ifft2(omegahat.*ky));
for i=1:nplots
    chg=1;
    % save old values
    uold=u; vold=v;  omegaold=omega; omegacheck=omega;
    omegahatold=omegahat; nonlinhatold=nonlinhat;
    while chg>tol
        %nonlinear {n+1,k}
        nonlinhat=fft2(u.*ifft2(omegahat.*kx)+v.*ifft2(omegahat.*ky));
       
        %Crank–Nicolson timestepping 
        omegahat=((1/dt + 0.5*(1/Re)*(kxx+kyy)).*omegahatold...
            -.5*(nonlinhatold+nonlinhat))...
            ./(1/dt -0.5*(1/Re)*(kxx+kyy));
        
        %compute \hat{\psi}^{n+1,k+1}    
        psihat=-omegahat./(kxx+kyy);  
        
        %NOTE: kxx+kyy has to be zero at the following points to avoid a
        % discontinuity. However, we suppose that the streamfunction has
        % mean value zero, so we set them equal to zero
        psihat(1,1)=0;                  
        psihat(Nx/2+1,Ny/2+1)=0;
        psihat(Nx/2+1,1)=0;
        psihat(1,Ny/2+1)=0;
        
        %computes {\psi}_x by differentiation via FFT
        dpsix = real(ifft2(psihat.*kx));  
        %computes {\psi}_y by differentiation via FFT
        dpsiy = real(ifft2(psihat.*ky));
        
        u=dpsiy;    %u^{n+1,k+1}
        v=-dpsix;   %v^{n+1,k+1}
        
        %\omega^{n+1,k+1}
        omega=ifft2(omegahat);
        % check for convergence
        chg=max(max(abs(omega-omegacheck))) 
        % store omega to check for convergence of next iteration
        omegacheck=omega;
    end     
     t(i+1)=t(i)+dt;  
     uexact_y=-2*pi*sin(2*pi*xx).*sin(2*pi*yy).*exp(-8*pi^2*t(i+1)/Re);
     vexact_x=2*pi*sin(2*pi*xx).*sin(2*pi*yy).*exp(-8*pi^2*t(i+1)/Re);
     omegaexact=vexact_x-uexact_y;
     figure(1); pcolor(omega);  xlabel x; ylabel y; 
     title Numerical; colorbar; drawnow;
     figure(2); pcolor(omegaexact); xlabel x; ylabel y; 
     title Exact; colorbar; drawnow;
     figure(3); pcolor(omega-omegaexact); xlabel x; ylabel y; 
     title Error; colorbar; drawnow;
end


