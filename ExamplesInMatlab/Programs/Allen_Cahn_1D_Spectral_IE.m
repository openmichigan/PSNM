%Solving 1D Allen-Cahn Eq using pseudo-spectral and Implicit/Explicit method
%u_t=u_{xx} + u - u^3
%where u-u^3 is treated explicitly and u_{xx} is treated implicitly
%BC = u(0)=0, u(2*pi)=0 (Periodic)
%IC=.25*sin(x); 
clear all; clc;

%Grid and Initial Data
N = 8000; h = 2*pi/N; x = h*(1:N); t = 0; 

dt = .001; %timestep size
epsilon= .001;

%initial conditions
v = .25*sin(x);   

%(ik) and (ik)^2 vectors
k=(1i*[0:N/2-1 0 -N/2+1:-1]);
k2=k.^2;

%setting up plot
tmax = 5; tplot = .2;
plotgap= round(tplot/dt); 
nplots = round(tmax/tplot);
data = [v; zeros(nplots,N)]; tdata = t;
  
for i = 1:nplots
    for n = 1:plotgap
        v_hat = fft(v); %converts to Fourier space
        vv = v.^3;      %computes nonlinear term in real space
        vv = fft(vv);   %converts nonlinear term to Fourier space  
        v_hat = (v_hat*(1/dt+1) - vv)./(1/dt-k2*epsilon); %Implicit/Explicit
        v = ifft(v_hat); %Solution back to real space
    end
    data(i+1,:) = real(v); %Records data each "plotgap"
    t=t+plotgap*dt; %Real time
    tdata = [tdata; t];
end

mesh(x,tdata,data), grid on, axis([-1 2*pi 0 tmax -1 1]),
view(-60,55), xlabel x, ylabel t, zlabel u