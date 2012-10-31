% A program to approximate an integral using the Monte Carlos method

% This program can be made much faster by using Matlab's matrix and vector
% operations, however to allow easy translation to other languages we have
% made it as simple as possible.

Numpoints=65536;   % number of random points

I2d=0; % Initialize value
I2dsquare=0; % initial variance
for n=1:Numpoints
    % generate random number drawn from a uniform distribution on (0,1)
    x=rand(1); 
    y=rand(1)*4;
    I2d=I2d+x^2+2*y^2;
    I2dsquare=I2dsquare+(x^2+2*y^2)^2;
end
% we sclae the integral by the total area and divide by the number of
% points used
I2d=I2d*4/Numpoints
% we also output an estimated error
I2dsquare=I2dsquare*4/Numpoints;
EstimError=4*sqrt( (I2d^2-I2dsquare)/Numpoints)