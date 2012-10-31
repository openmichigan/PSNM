%Solving 1D Allen-Cahn Eq using pseudo-spectral with Implicit Midpoint Rule
%u_t=u_{xx} + u - u^3
%where u-u^3 is treated explicitly and u_{xx} is treated implicitly
%BC = u(0)=0, u(2*pi)=0 (Periodic)
%IC=.25*sin(x); 
clear all; clc;

%Grid and Initial Data
N = 128; h = 2*pi/N; x = h*(1:N); t = 0; 

dt = .001;
epsilon= .001;

% initial conditions 
v=.25*sin(x);

%(ik) and (ik)^2 vectors
k=(1i*[0:N/2-1 0 -N/2+1:-1]);
k2=k.^2;

tol =10^-13; %tolerance

%setting up plot
tmax = 5; tplot = .2; clf, drawnow
plotgap= round(tplot/dt);
nplots = round(tmax/tplot);
data = [v; zeros(nplots,N)]; tdata = t;

  
for i = 1:nplots      
     for n=1:plotgap                   
        err=ones(1,N);
        vold = fft(v);
        vv = v.^3;
        vvold=fft(vv);
        while max(err)>tol %IMR iterations until tolerance is reached
            voldk=fft(v);
            vvk = v.^3;
            vvoldk=fft(vvk);
            vnewk= (vold.*(1/dt+k2*epsilon/2)...
                +(voldk-vvoldk)/2+(vold-vvold)/2)./(1/dt-k2*epsilon/2);
            err=abs(vnewk-voldk);
            v = ifft(vnewk);
        end              
     end
    data(i+1,:) = real(v); %records data
    t=t+plotgap*dt;
    tdata = [tdata; t];
end

%Plot
mesh(x,tdata,data), grid on,
view(-60,55), xlabel x, ylabel t, zlabel u;

