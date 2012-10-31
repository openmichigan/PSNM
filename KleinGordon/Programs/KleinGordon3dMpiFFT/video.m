% A program to create a video of the computed results

clear all; format compact, format short,
set(0,'defaultaxesfontsize',14,'defaultaxeslinewidth',.7,...
    'defaultlinelinewidth',2,'defaultpatchlinewidth',3.5);

% Load data
% Get coordinates
X=load('./xcoord.dat');
Y=load('./ycoord.dat');
Z=load('./zcoord.dat');
TIME=load('./tdata.dat');
En=load('./en.dat');
EnStr=load('./enstr.dat');
EnPot=load('./enpot.dat');
EnKin=load('./enkin.dat');

% find number of grid points
Nx=length(X);
Ny=length(Y);
Nz=length(Z);

% reshape coordinates to allow easy plotting
[xx,yy,zz]=meshgrid(X,Y,Z);


figure(5); clf; 
semilogy(TIME,En,'r',TIME,EnKin,'b:',TIME,EnPot,'g.',TIME,EnStr,'y+');
xlabel time; ylabel Energy; 
legend('Total','Kinetic','Potential','Strain');
saveas(5,'./EnerPlot.jpg','jpg');

% reshape coordinates to allow easy plotting
[xx,yy,zz]=meshgrid(X,Y,Z);

nplots=length(TIME);

for i =1:nplots
    %
    % Open file and dataset using the default properties.
    %
    FILE=['./data/u',num2str(10000000+i),'.datbin'];
    FILEPIC=['./data/pic',num2str(10000000+i),'.jpg'];
    FILEPICB=['./data/picB',num2str(10000000+i),'.jpg'];
    fid=fopen(FILE,'r');
    [fname,mode,mformat]=fopen(fid);
    u=fread(fid,Nx*Ny*Nz,'real*8');
    u=reshape(u,Nx,Ny,Nz);
    % close files
    fclose(fid);
    %
    % Plot data on the screen.
    %
    figure(2);clf; % coordinate slice to show plots on
	sx=[0]; sy=[0]; sz=[0];
	slice(xx,yy,zz,u,sx,sy,sz); colormap jet;
	title(['Time ',num2str(TIME(i))]); colorbar('location','EastOutside'); 
	drawnow; colorbar; frame=getframe(2); saveas(2,FILEPIC,'jpg');
	figure(4); clf; UP = u;
	p1 = patch(isosurface(X,Y,Z,UP,.1),...
    'FaceColor','yellow','EdgeColor','none');
	p2 = patch(isocaps(X,Y,Z,UP,.1),...
    'FaceColor','interp','EdgeColor','none');
	isonormals(UP,p1); lighting phong;
	xlabel('x'); ylabel('y'); zlabel('z'); 
	axis equal; axis square;  view(3); 
	title(['Time ',num2str(TIME(i))]); drawnow;
	saveas(4,FILEPICB,'jpg');
	
end

