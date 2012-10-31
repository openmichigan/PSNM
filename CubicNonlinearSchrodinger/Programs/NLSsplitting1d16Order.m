% A program to solve the nonlinear Schr\"{o}dinger equation using a
% splitting method. The numerical solution is compared to an exact
% solution.
% S. Blanes, F. Casas, P. Chartier and A. Murua
% "Optimized high-order splitting methods for some classes of parabolic
% equations"
% ArXiv pre-print 1102.1622v2
% Forthcoming Mathematics of Computation

clear all; format compact; format short;
set(0,'defaultaxesfontsize',30,'defaultaxeslinewidth',.7,...
    'defaultlinelinewidth',6,'defaultpatchlinewidth',3.7,...
    'defaultaxesfontweight','bold')

% set up grid
Lx = 20;        % period  2*pi * L
Nx = 16384;     % number of harmonics
Nt = 2000;      % number of time slices
dt = 0.25*pi/Nt;% time step
U=zeros(Nx,Nt/10);
method=3; % splitting method: 1 Strang, 2 CCDV10, 3 Blanes et al 2012

% initialise variables
x = (2*pi/Nx)*(-Nx/2:Nx/2 -1)'*Lx;          % x coordinate
kx = 1i*[0:Nx/2-1 0 -Nx/2+1:-1]'/Lx;        % wave vector

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

if method==1,
    %
    % Strang-Splitting
    %
    s=2;
    a=[1;0];
    b=[1/2;1/2];
    %
elseif method==2,
    %
    % Method of Castella, Chartier, Descombes and Vilmart 
    % BIT Numerical Analysis vol 49 pp 487-508, 2009
    %
    s=5;
    a=[1/4;1/4;1/4;1/4;0];
    b=[1/10-1i/30;4/15+2*1i/15;4/15-1i/5;4/15+2*1i/15;1/10-1i/30];
    %
elseif method==3,
    %
    % Method of Blanes, Casas, Chartier and Murua 2012
    %
    s=17;
    a=1/16*[1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;0];
    b=[0.028920177910074098791 - 0.005936580835725746103*1i;
       0.056654351383649876160 + 0.020841963949772627119*1i;
       0.067258385822722143569 - 0.039386393748812362460*1i;
       0.070333980553260772061 + 0.058952097930307840316*1i;
       0.077095100838099173580 - 0.038247636602014810025*1i;
       0.042022140317231098258 - 0.033116379859951038579*1i;
       0.050147397749937784280 + 0.061283684958324249562*1i;
       0.047750191909146143447 - 0.032332468814362628262*1i;
       0.119636547031757819706 + 0.015883426044923736862*1i;
       0.047750191909146143447 - 0.032332468814362628262*1i;
       0.050147397749937784280 + 0.061283684958324249562*1i;
       0.042022140317231098258 - 0.033116379859951038579*1i;
       0.077095100838099173580 - 0.038247636602014810025*1i;
       0.070333980553260772061 + 0.058952097930307840316*1i;
       0.067258385822722143569 - 0.039386393748812362460*1i;
       0.056654351383649876160 + 0.020841963949772627119*1i;
       0.028920177910074098791 - 0.005936580835725746103*1i];
end;


% solve pde and plot results
for n =2:Nt+1 
    for m=1:(s-1)     
        vna=exp(b(m)*1i*dt*kx.*kx).*v;
        una=ifft(vna);
        pot=(2*una.*conj(una));
        unb=exp(-1i*a(m)*(-1)*dt*pot).*una;
        v=fft(unb);      
    end
    v=exp(b(s)*1i*dt*kx.*kx).*v;
    u=ifft(v);
    t=(n-1)*dt;    
    if (mod(n,10)==0)
        tdata(n/10)=t;
        u=ifft(v);
        U(:,n/10)=u;
        uexact=...
            4*exp(1i*t)*(cosh(3*x)+3*exp(8*1i*t)*cosh(x))...
            ./(cosh(4*x)+4*cosh(2*x)+3*cos(8*t));
        figure(1); clf; plot(x,abs(u).^2); ...
            xlim([-0.5,0.5]); title(num2str(t));
        figure(2); clf; loglog(abs(v(1:Nx/2))); ...
            title('Fourier Coefficients');
        figure(3); clf; plot(x,abs(u-uexact).^2); ...
            xlim([-0.5,0.5]); title('error');
        drawnow;
        ma = fft(abs(u).^2);
        ma = ma(1);
        test = log10(abs(1-ma/ma0))
    end
end
figure(4); clf; mesh(tdata(1:(n-1)/10),x,abs(U(:,1:(n-1)/10)).^2);
xlim([0,t]);