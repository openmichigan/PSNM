%Solving 2D Allen-Cahn Eq using pseudo-spectral with Crank-Nicolson
%u_t= u_{xx}+u_{yy} + u - u^3
%BC = Periodic
%IC=v=sin(2*pi*x)+0.001*cos(16*pi*x;
clear all;clc;

%Grid
N = 64; h = 1/N; x = h*(1:N); 
y=x';
[xx,yy]=meshgrid(x,y);

dt = .01; 

%Initial Conditions
v=sin(2*pi*xx)+0.001*cos(16*pi*xx);
epsilon=.01;

%(ik) and (ik)^2 vectors in x and y direction
kx=(1i*[0:N/2-1 0 -N/2+1:-1]);
ky=(1i*[0:N/2-1 0 -N/2+1:-1]');
k2x=kx.^2;
k2y=ky.^2;

tol = 10^-8; %tolerance

[kxx,kyy]=meshgrid(k2x,k2y);
        
for n = 1:1000
    v_nl=v.^3; %calculates nonlinear term in real space
    %FFT for nonlinear and linear terms
    v_nl = fft2(v_nl);
    v_hat=fft2(v);
    err=1;
    while max(err)>tol %fixed point iterations until tolerance is reached
        v_nl_oldk=v.^3;  %oldk nonlinear term
        vvoldk = fft2(v_nl_oldk);
        voldk=fft2(v);       
        vnewk = (v_hat.*(1/dt+(uxx+uyy)*epsilon/2)...
                    +(voldk-v_nl_oldk)/2+(v_hat-v_nl)/2)...
              ./(1/dt-epsilon*(uxx/2+uyy/2)); %IMR Timestepping                
        err=sum(sum(abs(vnewk-voldk))); %Max error
        %Back to real space
        v=ifft2(vnewk);
    end
   %plot real part each timestep
   v=real(v);
   surf(v); title(num2str(n)); colormap hsv;  view(43,22); drawnow;
     
end
