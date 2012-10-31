%Solving Heat Equation using pseudospectral methods with Backwards Euler:
%u_t= \alpha*u_xx
%BC = u(0)=0 and u(2*pi)=0 (Periodic)
%IC=sin(x)
clear all; clc;

%Grid
N = 64; h = 2*pi/N; x = h*(1:N);     

% Initial conditions
v = sin(x);
alpha = .5; 
t = 0;          
dt = .001; %Timestep size

%(ik)^2 Vector
k=(1i*[0:N/2-1 0 -N/2+1:-1]);
k2=k.^2;

%Setting up Plot
tmax = 5; tplot = .1;
plotgap= round(tplot/dt); 
nplots = round(tmax/tplot);
data = [v; zeros(nplots,N)]; tdata = t;

  
for i = 1:nplots
    v_hat = fft(v); %Converts to fourier space
    for n = 1:plotgap
        v_hat = v_hat./(1-dt*alpha*k2); %Backwards Euler timestepping
    end
    v = ifft(v_hat); %Converts back to real Space
    data(i+1,:) = real(v); %Records data 
    t=t+plotgap*dt; %Records time
    tdata = [tdata; t];
end

%Plot using mesh
mesh(x,tdata,data), grid on, %axis([-1 2*pi 0 tmax -1 1]),
view(-60,55), xlabel x, ylabel t, zlabel u, zlabel u,
