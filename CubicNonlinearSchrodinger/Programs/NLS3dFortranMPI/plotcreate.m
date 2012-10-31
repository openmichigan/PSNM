% A program to plot the computed results

clear all; format compact, format short,
set(0,'defaultaxesfontsize',18,'defaultaxeslinewidth',.9,...
    'defaultlinelinewidth',3.5,'defaultpatchlinewidth',5.5);

% Load data
tdata=load('./tdata.dat');
x=load('./xcoord.dat');
y=load('./ycoord.dat');
z=load('./zcoord.dat');
Tsteps = length(tdata);

Nx = length(x); Nt = length(tdata);
Ny = length(y); Nz = length(z); 
fid=fopen('./ufinal.datbin','r');
[fname,mode,mformat]=fopen(fid);
u=fread(fid,Nx*Ny*Nz,'double',mformat);
u = reshape(u,Nx,Ny,Nz);

% Plot data
figure (1); clf ; UP = abs(u).^2; 
p1 = patch(isosurface(x,y,z,UP,.0025) ,...
	'FaceColor','yellow','EdgeColor','none'); 
p2 = patch(isocaps(x,y,z,UP,.0025) ,...
	'FaceColor','interp','EdgeColor','none'); 
isonormals(UP,p1); lighting phong;
xlabel('x'); ylabel('y'); zlabel('z');
axis equal; axis square; view(3); drawnow;
