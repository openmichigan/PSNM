%Solving Heat Equation using pseudo-spectral and Forward Euler
%u_t= \alpha*u_xx
%BC= u(0)=0, u(2*pi)=0
%IC=sin(x)
clear all; clc;

%Grid
N = 64;         %Number of steps
h = 2*pi/N;     %step size
x = h*(1:N);    %discretize x-direction

alpha = .5;      %Thermal Diffusivity constant
t = 0;          
dt = .001;

%Initial conditions 
v = sin(x);            
k=(1i*[0:N/2-1 0 -N/2+1:-1]);
k2=k.^2;

%Setting up Plot
tmax = 5; tplot = .1;
plotgap= round(tplot/dt); 
nplots = round(tmax/tplot);
data = [v; zeros(nplots,N)]; tdata = t;

  
for i = 1:nplots
    v_hat = fft(v);  %Fourier Space
    for n = 1:plotgap
        v_hat = v_hat+dt*alpha*k2.*v_hat; %FE timestepping
    end
    v = real(ifft(v_hat)); %Back to real space
    data(i+1,:) = v;
    t=t+plotgap*dt;
    tdata = [tdata; t]; %Time vector
end

%Plot using mesh
mesh(x,tdata,data), grid on,
view(-60,55), xlabel x, ylabel t, zlabel u, zlabel u