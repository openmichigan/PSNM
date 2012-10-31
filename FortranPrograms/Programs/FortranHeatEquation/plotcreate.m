% A Matlab program to plot the computed results

clear all; format compact, format short,
set(0,'defaultaxesfontsize',18,'defaultaxeslinewidth',.9,...
    'defaultlinelinewidth',3.5,'defaultpatchlinewidth',5.5);

% Load data
load('./u.dat');
load('./tdata.dat');
load('./xcoord.dat');
Tsteps = length(tdata);

Nx = length(xcoord); Nt = length(tdata);
 
u = reshape(u,Nx,Nt);

% Plot data
figure(3); clf; mesh(tdata,xcoord,u); xlabel t; ylabel x; zlabel('u');
