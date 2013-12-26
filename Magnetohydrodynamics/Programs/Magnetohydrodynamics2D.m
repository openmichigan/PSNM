% Numerical solution of the 2D incompressible Magnetohydrodynamics equations
% on a Square Domain [0,1]x[0,1] using a Fourier pseudo-spectral method
% and Implicit midpoint rule timestepping. A passive scalar is also advected 
% by the flow
%
%Periodic free-slip boundary conditions and Initial conditions:
    %u(x,y,0)=sin(2*pi*x)cos(2*pi*y)
    %v(x,y,0)=-cos(2*pi*x)sin(2*pi*y)
%Analytical Solution:
    %u(x,y,t)=sin(2*pi*x)cos(2*pi*y)exp(-8*pi^2*t/Re)
    %v(x,y,t)=-cos(2*pi*x)sin(2*pi*y)exp(-8*pi^2*t/Re)
% also consider an advection field
% \theta_t+u\theta_x+v\theta_y=Dtheta(\theta_{xx}+\theta_{yy})
clear all; format compact; format short; clc; clf;
 
Re=100;%Reynolds number
Rem=100; % Magnetic Reynolds number
Dtheta=0.01;%scalar diffusion constant
%grid
Nx=64; h=1/Nx; x=h*(1:Nx);
Ny=64; h=1/Ny; y=h*(1:Ny)';
[xx,yy]=meshgrid(x,y);
 
%initial conditions for velocity field
u=sin(2*pi*xx).*cos(2*pi*yy);
v=-cos(2*pi*xx).*sin(2*pi*yy);
u_y=-2*pi*sin(2*pi*xx).*sin(2*pi*yy);
v_x=2*pi*sin(2*pi*xx).*sin(2*pi*yy);
omega=v_x-u_y;
% initial magnetic potential field
alpha=sin(4*pi*xx).*cos(6*pi*yy);
% initial passive scalar field
theta=exp(-2.0*((4*pi*(xx-0.25)).^2+(4*pi*(yy-0.25)).^2));
dt=0.0025; t(1)=0; tmax=1.0;
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
%compute \hat{Phi}^{n+1,k+1}  
alphahat=fft2(alpha);
phihat=-alphahat./(kxx+kyy);  
 
%NOTE: kxx+kyy has to be zero at the following points to avoid a
% discontinuity. However, we suppose that the streamfunction has
% mean value zero, so we set them equal to zero
phihat(1,1)=0;                  
phihat(Nx/2+1,Ny/2+1)=0;
phihat(Nx/2+1,1)=0;
phihat(1,Ny/2+1)=0;
 
%computes {\psi}_x by differentiation via FFT
dphix = real(ifft2(phihat.*kx));  
%computes {\psi}_y by differentiation via FFT
dphiy = real(ifft2(phihat.*ky));
% components of magnetic field
bx=dphiy;    
by=-dphix;   
 
%compute \hat{\omega}^{n+1,k}
omegahat=fft2(omega); alphahat=fft2(alpha);
thetatotal(1)=sum(sum(theta.^2))*h*h;
for i=1:nplots
    chg=1;
    % save old values
    uold=u; vold=v;  omegaold=omega; omegacheck=omega;
    omegahatold=omegahat; thetaold=theta; thetacheck=theta;
    alphahatold=alphahat; alphaold=alpha; alphacheck=alpha;
    bxold=bx; byold=by;
    while chg>tol
        % Fluid field
        %nonlinear {n+1,k}
        nonlinhat=0.25*fft2((u+uold).*ifft2((omegahat+omegahatold).*kx)+...
                            (v+vold).*ifft2((omegahat+omegahatold).*ky)-...
                            (bx+bxold).*ifft2((alphahat+alphahatold).*kx)-...
                            (by+byold).*ifft2((alphahat+alphahatold).*ky));
 
        %Implicit midpoint rule timestepping 
        omegahat=((1/dt + 0.5*(1/Re)*(kxx+kyy)).*omegahatold-nonlinhat)...
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
 
        % magnetic field
        nonlinhat=0.25*fft2((u+uold).*ifft2((alphahat+alphahatold).*kx)+...
                            (v+vold).*ifft2((alphahat+alphahatold).*ky)-...
                            (bx+bxold).*ifft2((omegahat+omegahatold).*kx)-...
                            (by+byold).*ifft2((omegahat+omegahatold).*ky));
 
        %Implicit midpoint rule timestepping 
        alphahat=((1/dt + 0.5*(1/Re)*(kxx+kyy)).*alphahatold-nonlinhat)...
            ./(1/dt -0.5*(1/Rem)*(kxx+kyy));
 
        %compute \hat{\psi}^{n+1,k+1}    
        phihat=-alphahat./(kxx+kyy);  
 
        %NOTE: kxx+kyy has to be zero at the following points to avoid a
        % discontinuity. However, we suppose that the streamfunction has
        % mean value zero, so we set them equal to zero
        phihat(1,1)=0;                  
        phihat(Nx/2+1,Ny/2+1)=0;
        phihat(Nx/2+1,1)=0;
        phihat(1,Ny/2+1)=0;
 
        %computes {\psi}_x by differentiation via FFT
        dphix = real(ifft2(phihat.*kx));  
        %computes {\psi}_y by differentiation via FFT
        dphiy = real(ifft2(phihat.*ky));
 
        bx=dphiy;    %u^{n+1,k+1}
        by=-dphix;   %v^{n+1,k+1}
 
 
        thetax=0.5*ifft2(kx.*fft2(theta+thetaold));
        thetay=0.5*ifft2(ky.*fft2(theta+thetaold));
        theta=ifft2((fft2(thetaold-dt*0.5*((uold+u).*thetax+(vold+v).*thetay))+...
            dt*0.5*Dtheta*(kxx+kyy).*fft2(thetaold))./(1-dt*0.5*Dtheta*(kxx+kyy)));
 
        %\omega^{n+1,k+1}
        omega=ifft2(omegahat);
        % check for convergence
        chg=max(max(abs(omega-omegacheck))) + max(max(abs(theta-thetacheck)))+...
            max(max(abs(alpha-alphacheck)))
 
        % store omega and theta to check for convergence of next iteration
        omegacheck=omega; thetacheck=theta; alphacheck=alpha;
    end     
     t(i+1)=t(i)+dt;  
     thetatotal(i+1)=sum(sum(theta.^2))*h*h;
     figure(1); 
     subplot(2,2,1);
     pcolor(xx,yy,omega);  shading interp; xlabel x; ylabel y; 
     title(['Fluid Vorticity, Time ',num2str(t(i+1))]); colorbar; drawnow;
     subplot(2,2,2);
     pcolor(xx,yy,alpha); shading interp;  xlabel x; ylabel y; 
     title(['Alpha, Time ',num2str(t(i+1))]); colorbar; drawnow;
     subplot(2,2,3);
     pcolor(xx,yy,real(ifft2(psihat))); shading interp; xlabel x; ylabel y; 
     title(['Psi, Time ',num2str(t(i+1))]); colorbar; drawnow;
     subplot(2,2,4);
     pcolor(xx,yy,real(theta)); shading interp;  xlabel x; ylabel y; 
     title(['Theta, Time ',num2str(t(i+1))]); colorbar; drawnow;
end
figure(2); clf; 
subplot(2,1,1); plot(t,thetatotal,'r-'); xlabel time; ylabel('Total theta');
subplot(2,1,2); semilogy(t,abs(1-thetatotal./thetatotal(1)),'r-'); 
                xlabel time; ylabel('Total theta Error');

