% A program to plot the computed results for the 2D NLS equation

clear all; format compact, format short,
set(0,'defaultaxesfontsize',18,'defaultaxeslinewidth',.9,...
    'defaultlinelinewidth',3.5,'defaultpatchlinewidth',5.5);

% Load data
load('./ufinal.dat');
load('./tdata.dat');
load('./ycoord.dat');
load('./xcoord.dat');

Ny = length(ycoord); Nx = length(xcoord); Nt = length(tdata);
 
ufinal = reshape(ufinal,Nx,Ny);

% Plot data
figure(3); clf; mesh(xcoord,ycoord,ufinal); xlabel x; ylabel y; zlabel('|u|^2');