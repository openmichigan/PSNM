% Please use an appropriate header since this file does not solve the 2nd
% order wave equation. Please also ensure that you credit any sources of
% code you have used/modified to get your program.

%p6.m - 2nd Order Wave equation
% Grid, variable coefficient, and initial data:
clear all;
clc;
N = 64; h = 1/N; x = h*(1:N); t = 0; dt = .001; %When dt=.01 it blows up for N=64
y=x';
[xx,yy]=meshgrid(x,y);

vv = exp((-200*(xx-.4).^2 + -100*(yy-.4).^2));
kx=(1i*[0:N/2-1 0 -N/2+1:-1]);
ky=(1i*[0:N/2-1 0 -N/2+1:-1]');
k2x=kx.^2;
k2y=ky.^2;
ii=1:N;
for n = 1:5000
    v_hat=fft2(vv);
    for m= 1:N
        vx = v_hat(m,:);
        uxx(m,:)= k2x(ii).*vx;
    end
    for j= 1:N
        vy = v_hat(:,j);
        uyy(:,j)= k2y(ii).*vy;
    end
    vnew = v_hat + dt*(uyy+uxx);
    vv=ifft2(vnew);   
    surf(vv); title(num2str(n)); axis([0 65 0 65 -.2 .5]); colormap hsv; drawnow;
     
end