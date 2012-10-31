% A program to approximate an integral

clear all; format compact; format short;

nx=1000;     % number of points in x
xend=1;     % last discretization point
xstart=0;   % first discretization point
dx=(xend-xstart)/(nx-1);    % size of each x sub-interval

ny=4000;     % number of points in y
yend=4;     % last discretization point
ystart=0;   % first discretization point
dy=(yend-ystart)/(ny-1);    % size of each y sub-interval

% create vectors with points for x and y
for i=1:nx
    x(i)=xstart+(i-1)*dx;
end
for j=1:ny
    y(j)=ystart+(j-1)*dy;
end

% Approximate the integral by a sum
I2d=0;
for i=1:nx
    for j=1:ny
        I2d=I2d+(x(i)^2+2*y(j)^2)*dy*dx;
    end
end
% print out final answer
I2d