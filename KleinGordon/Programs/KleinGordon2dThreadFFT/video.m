% A program to create a video of the computed results

clear all; format compact, format short,
set(0,'defaultaxesfontsize',14,'defaultaxeslinewidth',.7,...
    'defaultlinelinewidth',2,'defaultpatchlinewidth',3.5);

% Load data
% Get coordinates
X=load('./xcoord.dat');
Y=load('./ycoord.dat');
TIME=load('./tdata.dat');
% find number of grid points
Nx=length(X);
Ny=length(Y);

% reshape coordinates to allow easy plotting
[xx,yy]=ndgrid(X,Y);

nplots=length(TIME);

for i =1:nplots
    %
    % Open file and dataset using the default properties.
    %
    FILE=['./data/u',num2str(9999999+i),'.datbin'];
    FILEPIC=['./data/pic',num2str(9999999+i),'.jpg'];
    fid=fopen(FILE,'r');
    [fname,mode,mformat]=fopen(fid);
    u=fread(fid,Nx*Ny,'real*8');
    u=reshape(u,Nx,Ny);
    % close files
    fclose(fid);
    %
    % Plot data on the screen.
    %
    figure(2);clf; mesh(xx,yy,real(u)); xlabel x; ylabel y;
    title(['Time ',num2str(TIME(i))]); colorbar; axis square; 
    drawnow; frame=getframe(2); saveas(2,FILEPIC,'jpg');
end

