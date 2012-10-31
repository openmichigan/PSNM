%2D Heat Equation
%u_t=\alpha(u_{xx}+u_{yy})
%BC= Periodic in x and y directions
%IC= exp((-200*(xx-.5).^2 + -100*(yy-.5).^2))
clear all; clc;

% Grid and Initial Data
N = 64; h = 1/N; x = h*(1:N); t = 0; dt = .001;

%x and y meshgrid
y=x';
[xx,yy]=meshgrid(x,y);

%Initial Conditions
vv = exp((-200*(xx-.5).^2 + -100*(yy-.5).^2));

%(ik)^2 vectors in x and y directions
kx=(1i*[0:N/2-1 0 -N/2+1:-1]);
ky=(1i*[0:N/2-1 0 -N/2+1:-1]');
k2x=kx.^2;
k2y=ky.^2;

[kxx,kyy]=meshgrid(k2x,k2y);
       
for n = 1:1000
    v_hat=fft2(vv);        %FFT 
    vnew = v_hat./(1-dt*(kyy+kxx));      %Backwards Euler
    vv=ifft2(vnew);        %Back to real space
    surf(vv); title(num2str(n)); axis([0 N+1 0 N+1 -.2 .5]); %Plot in real 
    colormap hsv; drawnow; 
end