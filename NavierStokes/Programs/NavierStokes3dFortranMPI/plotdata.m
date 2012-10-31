% A program to create a plot of the computed results
% from the 3D Fortran Navier-Stokes solver
clear all; format compact; format short;
set(0,'defaultaxesfontsize',30,'defaultaxeslinewidth',.7,...
    'defaultlinelinewidth',6,'defaultpatchlinewidth',3.7,...
    'defaultaxesfontweight','bold')
% Load data
% Get coordinates
tdata=load('./data/tdata.dat');
x=load('./data/xcoord.dat');
y=load('./data/ycoord.dat');
z=load('./data/zcoord.dat');
nplots = length(tdata);

Nx = length(x); Nt = length(tdata);
Ny = length(y); Nz = length(z); 

% reshape coordinates to allow easy plotting
[xx,yy,zz]=meshgrid(x,y,z);

for i =1:nplots
    %
    % Open file and dataset using the default properties.
    %
    FILEX=['./data/omegax',num2str(9999999+i),'.datbin'];
    FILEY=['./data/omegay',num2str(9999999+i),'.datbin'];
    FILEZ=['./data/omegaz',num2str(9999999+i),'.datbin'];
    FILEPIC=['./data/pic',num2str(9999999+i),'.jpg'];
    fid=fopen(FILEX,'r');
    [fname,mode,mformat]=fopen(fid);
    omegax=fread(fid,Nx*Ny*Nz,'real*8');
    omegax=reshape(omegax,Nx,Ny,Nz);
    fid=fopen(FILEY,'r');
    [fname,mode,mformat]=fopen(fid);
    omegay=fread(fid,Nx*Ny*Nz,'real*8');
    omegay=reshape(omegay,Nx,Ny,Nz);
    fid=fopen(FILEZ,'r');
    [fname,mode,mformat]=fopen(fid);
    omegaz=fread(fid,Nx*Ny*Nz,'real*8');
    omegaz=reshape(omegaz,Nx,Ny,Nz);
    % close files
    fclose(fid);
    %
    % Plot data on the screen.
    %
    omegatot=omegax.^2+omegay.^2+omegaz.^2;
    figure(100); clf; 
    subplot(2,2,1); title(['omega x ',num2str(tdata(i))]);
    p1 = patch(isosurface(x,y,z,omegax,.0025),...
            'FaceColor','interp','EdgeColor','none','FaceAlpha',0.3);
    p2 = patch(isocaps(x,y,z,omegax,.0025),...
            'FaceColor','interp','EdgeColor','none','FaceAlpha',0.1);
        isonormals(omegax,p1); lighting phong;
    xlabel('x'); ylabel('y'); zlabel('z'); 
    axis equal; axis square; view(3); colorbar;
    subplot(2,2,2); title(['omega y ',num2str(tdata(i))]);
    p1 = patch(isosurface(x,y,z,omegay,.0025),...
            'FaceColor','interp','EdgeColor','none','FaceAlpha',0.3);
    p2 = patch(isocaps(x,y,z,omegay,.0025),...
            'FaceColor','interp','EdgeColor','none','FaceAlpha',0.1);
        isonormals(omegay,p1); lighting phong;
    xlabel('x'); ylabel('y'); zlabel('z'); 
    axis equal; axis square; view(3); colorbar;
    subplot(2,2,3); title(['omega z ',num2str(tdata(i))]);
    p1 = patch(isosurface(x,y,z,omegaz,.0025),...
            'FaceColor','interp','EdgeColor','none','FaceAlpha',0.3);
    p2 = patch(isocaps(x,y,z,omegaz,.0025),...
            'FaceColor','interp','EdgeColor','none','FaceAlpha',0.1);
        isonormals(omegaz,p1); lighting phong;
    xlabel('x'); ylabel('y'); zlabel('z'); 
    axis equal; axis square; view(3); colorbar;
    subplot(2,2,4); title(['|omega|^2 ',num2str(tdata(i))]);
    p1 = patch(isosurface(x,y,z,omegatot,.0025),...
            'FaceColor','interp','EdgeColor','none','FaceAlpha',0.3);
    p2 = patch(isocaps(x,y,z,omegatot,.0025),...
            'FaceColor','interp','EdgeColor','none','FaceAlpha',0.1);
        isonormals(omegatot,p1); lighting phong;
    xlabel('x'); ylabel('y'); zlabel('z'); colorbar;
    axis equal; axis square; view(3);
	saveas(100,FILEPIC);
	
end

