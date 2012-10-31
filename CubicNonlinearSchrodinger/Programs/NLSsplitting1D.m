% A program to solve the nonlinear Schr\"{o}dinger equation using a
% splitting method

clear all; format compact; format short;
set(0,'defaultaxesfontsize',30,'defaultaxeslinewidth',.7,...
    'defaultlinelinewidth',6,'defaultpatchlinewidth',3.7,...
    'defaultaxesfontweight','bold')

Lx = 20;            % period  2*pi * L
Nx = 16384;         % number of harmonics
Nt = 1000;          % number of time slices
dt = 0.25*pi/Nt;    % time step
U=zeros(Nx,Nt/10);

Es = -1; % focusing or defocusing parameter

% initialise variables
x = (2*pi/Nx)*(-Nx/2:Nx/2 -1)'*Lx;      % x coordinate
kx = 1i*[0:Nx/2-1 0 -Nx/2+1:-1]'/Lx;    % wave vector
k2x = kx.^2;                       % square of wave vector
% initial conditions
t=0; tdata(1)=t;
u=4*exp(1i*t)*(cosh(3*x)+3*exp(8*1i*t)*cosh(x))...
    ./(cosh(4*x)+4*cosh(2*x)+3*cos(8*t));
v=fft(u);
figure(1); clf; plot(x,u);xlim([-2,2]); drawnow;
U(:,1)=u; 

% mass
ma = fft(abs(u).^2);
ma0 = ma(1);

% solve pde and plot results
for n =2:Nt+1
    
    vna=exp(0.5*1i*dt*k2x).*v;
    una=ifft(vna);
    pot=2*(una.*conj(una));
    unb=exp(-1i*Es*dt*pot).*una;
    vnb=fft(unb);
    v=exp(0.5*1i*dt*k2x).*vnb;
    t=(n-1)*dt;
    
    if (mod(n,10)==0)
        tdata(n/10)=t;
        u=ifft(v);
        U(:,n/10)=u;
        uexact=4*exp(1i*t)*(cosh(3*x)+3*exp(8*1i*t)*cosh(x))...
            ./(cosh(4*x)+4*cosh(2*x)+3*cos(8*t));
        figure(1); clf; plot(x,abs(u).^2); ...
            xlim([-0.5,0.5]); title(num2str(t));
        figure(2); clf; plot(x,abs(u-uexact).^2);...
            xlim([-0.5,0.5]); title(num2str(t));
        drawnow;
        ma = fft(abs(u).^2);
        ma = ma(1);
        test = log10(abs(1-ma/ma0))
    end
end
figure(3); clf; mesh(tdata(1:(n-1)/10),x,abs(U(:,1:(n-1)/10)).^2);