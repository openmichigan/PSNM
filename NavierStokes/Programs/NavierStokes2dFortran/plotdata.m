% A program to create a plot of the computed results
% from the 2D Matlab Navier-Stokes solver

clear all; format compact, format short,
set(0,'defaultaxesfontsize',14,'defaultaxeslinewidth',.7,...
    'defaultlinelinewidth',2,'defaultpatchlinewidth',3.5);

% Load data
% Get coordinates
X=load('xcoord.dat');
Y=load('ycoord.dat');
% find number of grid points
Nx=length(X);
Ny=length(Y);

% reshape coordinates to allow easy plotting
[xx,yy]=ndgrid(X,Y);

%
% Open file and dataset using the default properties.
%
FILENUM=['omegafinal.datbin'];
FILEEXA=['omegaexactfinal.datbin'];
fidnum=fopen(FILENUM,'r');
[fnamenum,modenum,mformatnum]=fopen(fidnum);
fidexa=fopen(FILEEXA,'r');
[fnameexa,modeexa,mformatexa]=fopen(fidexa);
Num=fread(fidnum,Nx*Ny,'double',mformatnum);
Exa=fread(fidexa,Nx*Ny,'double',mformatexa);
Num=reshape(Num,Nx,Ny);
Exa=reshape(Exa,Nx,Ny);
% close files
fclose(fidnum);
fclose(fidexa);
%
% Plot data on the screen.
%
figure(2);clf; 
subplot(3,1,1); contourf(xx,yy,Num);
title(['Numerical Solution ']); 
colorbar; axis square; 
subplot(3,1,2); contourf(xx,yy,Exa); 
title(['Exact Solution ']);
colorbar; axis square;
subplot(3,1,3); contourf(xx,yy,Exa-Num); 
title(['Error']);
colorbar; axis square;
drawnow; 