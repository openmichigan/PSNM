%Solving 2D Allen-Cahn Eq using pseudo-spectral with Implicit/Explicit
%u_t= epsilon(u_{xx}+u_{yy}) + u - u^3
%where u-u^3 is treated explicitly and epsilon(u_{xx} + u_{yy}) is treated implicitly
%BC = Periodic
%IC=v=sin(2*pi*x)+0.001*cos(16*pi*x;
clear all; clc;

%Grid
N = 256; h = 1/N; x = h*(1:N); 
dt = .01;

%x and y meshgrid
y=x';
[xx,yy]=meshgrid(x,y);

%initial conditions
v=sin(2*pi*xx)+0.001*cos(16*pi*xx);
epsilon=.01;

%(ik) and (ik)^2 vectors in x and y direction
kx=(1i*[0:N/2-1 0 -N/2+1:-1]);
ky=(1i*[0:N/2-1 0 -N/2+1:-1]');
k2x=kx.^2;
k2y=ky.^2;

[kxx,kyy]=meshgrid(k2x,k2y);
        
for n = 1:500   
    v_nl=v.^3;  %calculates nonlinear term in real space
    %FFT for linear and nonlinear term
    v_nl = fft2(v_nl);
    v_hat=fft2(v);
    vnew=(v_hat*(1+1/dt)-v_nl)./ ...
       (-(kxx+kyy)*epsilon+1/dt); %Implicit/Explicit timestepping
    %converts to real space in x-direction
    v=ifft2(vnew);   
    %Plots each timestep
    mesh(v); title(['Time ',num2str(n)]); axis([0 N 0 N -1 1]); 
    xlabel x; ylabel y; zlabel u;
    view(43,22); drawnow;  
end