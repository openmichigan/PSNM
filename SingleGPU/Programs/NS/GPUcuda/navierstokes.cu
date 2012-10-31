/*--------------------------------------------------------------------
*
* PURPOSE
*
* This program numerically solves the 2D incompressible Navier-Stokes
* on a Square Domain [0,1]x[0,1] using pseudo-spectral methods and
* Crank-Nicolson timestepping. The numerical solution is compared to
* the exact Taylor-Green Vortex Solution. 
*
* Periodic free-slip boundary conditions and Initial conditions:
*	u(x,y,0)=sin(2*pi*x)cos(2*pi*y)
*	v(x,y,0)=-cos(2*pi*x)sin(2*pi*y)
* Analytical Solution (subscript denote derivatives):
*	u(x,y,t)=sin(2*pi*x)cos(2*pi*y)exp(-8*pi^2*t/Re)
*	v(x,y,t)=-cos(2*pi*x)sin(2*pi*y)exp(-8*pi^2*t/Re)
*   u_y(x,y,t)=-2*pi*sin(2*pi*x)sin(2*pi*y)exp(-8*pi^2*t/Re)
*	v_x(x,y,t)=2*pi*sin(2*pi*x)sin(2*pi*y)exp(-8*pi^2*t/Re)
*	omega=v_x-u_y
*
* AUTHORS
*
* B. Cloutier, B.K. Muite, P. Rigge
* 4 June 2012
*
* .. Scalars ..
*  Nx				= number of modes in x - power of 2 for FFT
*  Ny				= number of modes in y - power of 2 for FFT
*  nplots			= number of plots produced
*  plotgap			= number of timesteps inbetween plots
*  Re 				= Reynold's number
*  dt				= timestep size 
*  tol				= determines when convergences is reached
*  i				= loop counter in x direction
*  j				= loop counter in y direction
*  n				= loop counter for timesteps between plots	
*  nn 				= loop counter for plots
*  chg				= error at each iteration	
*  max				= maximum error
*  pi				= value of pi
*  xsize    		= size of real arrays in x direction
*  ysize    		= size of real arrays in y direction
*  gridsize 		= size of array for complex data
*  start_time		= start time of computation
*  end_time			= end time of evaluation
*  pland2z			= Forward 2d fft plan  (CUFFT)
*  planz2d 			= Backward 2d fft plan (CUFFT)
*  nThreads			= Number of threads for GPU to use
*  nBlocksR			= number of blocks for GPU to use for real arrays
*  nBlocksC			= number of blocks for GPU to use for complex arrays
*
* .. Arrays on CPU ..
*
*  u				= velocity in x direction
*  v				= velocity in y direction
*  omeg				= vorticity	in real space
*  omegold			= vorticity in real space at previous
*						iterate
*  omegexact		= taylor-green vorticity at
*						at final step
*  x				= x locations
*  y				= y locations
*
* .. Arrays on GPU ..
*
*  u_d				= velocity in x direction
*  v_d				= velocity in y direction
*  omeg_d			= vorticity	in real space
*  omegold			= vorticity in real space at previous
*						iterate
*  omegoldhat_d		= 2D Fourier transform of vorticity at previous
*						iterate
*  nlhat_d			= nonlinear term in Fourier space
*  nloldhat_d		= nonlinear term in Fourier space
*						at previous iterate
*  omegexact		= taylor-green vorticity at
*						at final step
*  psihat_d			= 2D Fourier transform of streamfunction
*						at next iteration
*  omegcheck_d		= store of vorticity at previous iterate
*  temp1_d 			= temporary real space used for
*						calculations.
*  temp2_d 			= temporary complex space used for
*						calculations. 
*  kx_d				= fourier frequencies in x direction
*  ky_d				= fourier frequencies in y direction
*  x_d				= x locations
*  y_d				= y locations
*
* REFERENCES
*
* ACKNOWLEDGEMENTS
*
* The format for the complex data types and style has followed examples on 
* the Nvidia website
*
* ACCURACY
*		
* ERROR INDICATORS AND WARNINGS
*
* FURTHER COMMENTS
* Check that the initial iterate is consistent with the
* boundary conditions for the domain specified
*--------------------------------------------------------------------
* External libraries required
*       Cuda FFT
*/


#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>

#include <cufft.h>
#include <cutil_inline.h>

typedef double2 cmplx;
typedef double real;
typedef cufftDoubleReal cufftReal_;
typedef cufftDoubleComplex cufftCmplx;
#define cufftMyExecF cufftExecD2Z
#define cufftMyExecB cufftExecZ2D
#define CUFFT_TFORM_FORWARD  CUFFT_D2Z
#define CUFFT_TFORM_BACKWARD CUFFT_Z2D

static __device__ __host__ inline cmplx cmul(cmplx a, cmplx b){
	cmplx c;
	c.x = a.x * b.x - a.y * b.y;
	c.y = a.x * b.y + a.y * b.x;
	return c;
}

static __device__ inline cmplx cadd (cmplx a, cmplx b) {
  cmplx c;
  c.x = a.x + b.x;
  c.y = a.y + b.y;
  return c;
}

static __device__ inline cmplx csub (cmplx a, cmplx b) {
  cmplx c;
  c.x = a.x - b.x;
  c.y = a.y - b.y;
  return c;
}

static __device__ inline cmplx cscale(cmplx a, real b) {
	cmplx c;
	c.x = a.x * b;
	c.y = a.y * b;
	return c;
}

void initialdata (int Nx, int Ny, real pi, real *x, real *y, real *omeg) {
	int i, j;
	for(j=0; j<Ny; j++){
		for(i=0; i<Nx; i++){
			omeg[Nx*j + i] = 4.0*pi*sin(2.0*pi*x[i])*sin(2.0*pi*y[j]);
		}
	}
}

static __global__ void nonlin1(int Ny, real *kx_d, cmplx *omeghat_d, cmplx *temp1_d) {
  int idx = blockIdx.x * blockDim.x + threadIdx.x;
  int i;
  i = idx / (Ny/2+1);	
  cmplx c;
  c.x=(real)0.0;
  c.y=kx_d[i];
  temp1_d[idx]=cmul(omeghat_d[idx], c);
}

static __global__ void nonlin2(real *u_d, real *temp2_d, real *nl_d) {
	int idx = blockIdx.x * blockDim.x + threadIdx.x;
	nl_d[idx]=u_d[idx]*temp2_d[idx];
}

static __global__ void nonlin3(int Ny, real *ky_d, cmplx *omeghat_d, cmplx *temp1_d) {
	int idx = blockIdx.x * blockDim.x + threadIdx.x;
	int j;
	cmplx c;
	j = idx % (Ny/2+1);	
	c.x=(real)0.0;
	c.y=ky_d[j];
	temp1_d[idx]=cmul(omeghat_d[idx], c);
}

static __global__ void nonlin4(int Nx, int Ny, real *nl_d, real *v_d, real *temp2_d) {
	int idx = blockIdx.x * blockDim.x + threadIdx.x;
	real scale= (real) 1./(real)(Nx*Ny);
	nl_d[idx]=(nl_d[idx]+v_d[idx]*temp2_d[idx])*scale;
}

static __global__ void nextstep1(int Nx, int Ny, real dt, real Re, real *kx_d, real *ky_d, cmplx *omeghat_d, cmplx *omegoldhat_d, cmplx *nloldhat_d, cmplx *nlhat_d) {
	int idx = blockIdx.x * blockDim.x + threadIdx.x;
	int i,j;
	i = idx / (Ny/2+1);
	j = idx % (Ny/2+1);	
	
	real coef1 = (real) 1.*(-kx_d[i]*kx_d[i] - ky_d[j]*ky_d[j]);
	real coef2 = (real) 1./(dt);
	real coef3 = (real) 1./(Re);
	real coef4 = (real) 1.*(coef2 + (real)0.5*coef3*coef1);
	
	cmplx term1 = cscale (omegoldhat_d[idx],coef4);
	cmplx term2 = cadd (nloldhat_d[idx],nlhat_d[idx]);
	cmplx term3 = cscale (term2, (real)0.5);
	real coef5 = (real) 1./(coef2 - (real)0.5*coef3*coef1);
	omeghat_d[idx] = cscale(csub (term1,term3), coef5);

}

static __global__ void nextstep2(int Nx, int Ny, real *kx_d, real *ky_d, cmplx *omeghat_d, cmplx *psihat_d) {
	int idx = blockIdx.x * blockDim.x + threadIdx.x;
	int i,j;
	i = idx / (Ny/2+1);
	j = idx % (Ny/2+1);	

	real coef6 = (real) -1./(-kx_d[i]*kx_d[i] - ky_d[j]*ky_d[j] + pow((real)0.10,14));
	psihat_d[idx] = cscale(omeghat_d[idx],coef6);	
}

static __global__ void nextstep3(int Nx, int Ny, cmplx *psihat_d, real *kx_d, cmplx *temp1_d) {
	int idx = blockIdx.x * blockDim.x + threadIdx.x;
	int i;
	cmplx c;
	i = idx / (Ny/2+1);	
	c.x=(real)0.0;
	c.y=-kx_d[i]/(real)(Nx*Ny);
	temp1_d[idx]= cmul(c,psihat_d[idx]);
}
	
static __global__ void nextstep4(int Nx, int Ny, cmplx *psihat_d, real *ky_d, cmplx *temp1_d){
	int idx = blockIdx.x * blockDim.x + threadIdx.x;
	int j;
	cmplx c;
	j = idx % (Ny/2+1);
	c.x=(real)0.0;
	c.y=ky_d[j]/(real)(Nx*Ny);
	temp1_d[idx]=cmul(c,psihat_d[idx]);
}
			
static __global__ void copyRealArray(real *lhs, real *rhs) {
	int idx = blockIdx.x * blockDim.x + threadIdx.x;
	lhs[idx]=rhs[idx];
}

static __global__ void copyCmplxArray(cmplx *lhs, cmplx *rhs) {
	int idx = blockIdx.x * blockDim.x + threadIdx.x;
	lhs[idx].x=rhs[idx].x;
	lhs[idx].y=rhs[idx].y;
}

static __global__ void checkConvergence1(int Nx, int Ny, real *omeg_d, real *omegcheck_d) {
	int idx = blockIdx.x * blockDim.x + threadIdx.x;
	real scale=(real)1.0/(real)(Nx*Ny);
	omegcheck_d[idx] = ((omeg_d[idx]-omegcheck_d[idx])*scale);
	omegcheck_d[idx] = omegcheck_d[idx]*omegcheck_d[idx]; 
}


void savedata (int Nx, int Ny, int nplot, real *omeg) {
  FILE *f = NULL;
  char nameconfig1[128];
  nameconfig1[0]='\0';
  sprintf (nameconfig1, "data/omeg%.10d.datbin", nplot);
  f = fopen (nameconfig1, "wb");
  real *omegscale;
  omegscale = (real*)malloc (Nx * Ny * sizeof(real));
  int i;
  for (i=0; i<Nx*Ny; i++) {
    omegscale[i] = omeg[i]*1./(real)(Nx*Ny);
  }
  fwrite (omegscale, sizeof(real), Nx*Ny, f);
  fclose (f);
}

int main(){
	// declare variables
	cufftHandle planz2d, pland2z;	
	real chg, max;
	real pi;
	int Nx=1024;
	int Ny=1024;
	int nplots=1;
	int plotgap=20;
	real Re = 1.e0;
	real tol = 10.e-10;
	real dt = 0.000125;
	int nThreads;
	int nBlocksR;
	int nBlocksC;
	real *kx,*ky;
	real *x, *y, *omeg, *omegexact;
	int i, j, n, nn;
	struct timeval start_time, end_time;
	// declare variables for GPU
	real *u_d, *v_d;
	real *omegcheck_d, *omeg_d, *nl_d, *temp2_d;
	cmplx *omegoldhat_d, *nloldhat_d, *omeghat_d, *nlhat_d, *psihat_d, *temp1_d;
	real *kx_d, *ky_d;
	real *x_d, *y_d;
	size_t xsize	= Nx * sizeof(real);
	size_t ysize	= Ny * sizeof(real);
	size_t gridsize = Nx * (Ny/2+1)*sizeof(cmplx);

	pi = 3.14159265358979323846264338327950288419716939937510;
	nThreads= (min(512, Nx));
	nBlocksR= (Nx * Ny / nThreads);
	nBlocksC= (Nx * (Ny/2+1) / nThreads);
	
	printf("Program starting\n");
	printf("Grid: %d X %d\n",Nx,Ny);
	printf("dt: %lf\n",dt);
	
	kx=(real*) malloc(xsize);
	ky=(real*) malloc(ysize);
	x=(real*) malloc(xsize);
	y=(real*) malloc(ysize);
	omeg=(real*) malloc(Nx * Ny * sizeof(real));
	omegexact=(real*) malloc(Nx * Ny * sizeof(real));

	printf("Allocated CPU arrays\n");
	
	cutilSafeCall(cudaMalloc((void**)&u_d, Nx * Ny * sizeof(cufftReal_)));
	cutilSafeCall(cudaMalloc((void**)&v_d, Nx * Ny * sizeof(cufftReal_)));
	cutilSafeCall(cudaMalloc((void**)&omegcheck_d, Nx * Ny * sizeof(cufftReal_)));
	cutilSafeCall(cudaMalloc((void**)&omeg_d, Nx * Ny * sizeof(cufftReal_)));
	cutilSafeCall(cudaMalloc((void**)&nl_d, Nx * Ny * sizeof(cufftReal_)));
	cutilSafeCall(cudaMalloc((void**)&temp2_d, Nx * Ny * sizeof(cufftReal_)));
	cutilSafeCall(cudaMalloc((void**)&omegoldhat_d, gridsize));
	cutilSafeCall(cudaMalloc((void**)&nloldhat_d, gridsize));
	cutilSafeCall(cudaMalloc((void**)&omeghat_d, gridsize));
	cutilSafeCall(cudaMalloc((void**)&nlhat_d, gridsize));
	cutilSafeCall(cudaMalloc((void**)&psihat_d, gridsize));
	cutilSafeCall(cudaMalloc((void**)&temp1_d, gridsize));
	cutilSafeCall(cudaMalloc((void**)&kx_d, Nx * sizeof(cufftReal_)));
	cutilSafeCall(cudaMalloc((void**)&ky_d, Ny * sizeof(cufftReal_)));
	cutilSafeCall(cudaMalloc((void**)&x_d, Nx * sizeof(cufftReal_)));
	cutilSafeCall(cudaMalloc((void**)&y_d, Ny * sizeof(cufftReal_)));
	printf("Allocated GPU arrays\n");

	cufftSafeCall(cufftPlan2d(&pland2z, Nx, Ny, CUFFT_TFORM_FORWARD));
	cufftSafeCall(cufftPlan2d(&planz2d, Nx, Ny, CUFFT_TFORM_BACKWARD));
	printf("Setup FFTs\n");
	
	// setup fourier frequencies
	for(i=0; i<Nx/2; i++)
		kx[i]=2.0*pi*(real)i; 			
	kx[Nx/2]=0;	
	for(i=0; i<Nx/2-1; i++)
		kx[Nx/2+1+i] = -kx[Nx/2-1-i];
	for(i=0; i<Nx; i++)
		x[i]=(real)i/(real)Nx; 
		
	for(j=0; j<Ny/2; j++)
		ky[j]=2.0*pi*(real)j;  			
	ky[Ny/2]=0.0;
	for(j=0; j<Ny/2-1; j++)
		ky[Ny/2+1+j]=-ky[Ny/2-1-j];
	for(j=0; j<Ny; j++)
		y[j]=(real)j/(real)Ny; 
		
	cutilSafeCall(cudaMemcpy(kx_d,kx,Nx*sizeof(real),cudaMemcpyHostToDevice));
	cutilSafeCall(cudaMemcpy(x_d,x,Nx*sizeof(real),cudaMemcpyHostToDevice));
	cutilSafeCall(cudaMemcpy(ky_d,ky,Ny*sizeof(real),cudaMemcpyHostToDevice));
	cutilSafeCall(cudaMemcpy(y_d,y,Ny*sizeof(real),cudaMemcpyHostToDevice));
	
	printf("Setup grid and fourier frequencies\n");

	//!!!!!!!!!!!!!!
	//!initial data!
	//!!!!!!!!!!!!!!
	initialdata(Nx,Ny,pi,x,y,omeg);
	cutilSafeCall(cudaMemcpy(omeg_d,omeg,Ny*Nx*sizeof(real),cudaMemcpyHostToDevice));
	printf("Copied initial data to device\n");
	
	copyRealArray <<< nBlocksR, nThreads >>>(omegcheck_d, omeg_d);
	cutilCheckMsg("Kernel execution failed: [ copyRealArray ]");
	cufftSafeCall(cufftMyExecF(pland2z,(cufftReal_ *)omeg_d,(cufftCmplx *)omeghat_d));

	nextstep2 <<< nBlocksC, nThreads >>> (Nx, Ny, kx_d, ky_d, omeghat_d, psihat_d);
	cutilCheckMsg("Kernel execution failed: [ nextstep2 ]");
 				
	nextstep3 <<< nBlocksC, nThreads >>>(Nx,Ny,psihat_d,kx_d,temp1_d);
	cutilCheckMsg("Kernel execution failed: [ nextstep3 ]");
	cufftSafeCall(cufftMyExecB(planz2d,(cufftCmplx *)temp1_d,(cufftReal_ *)v_d));

	nextstep4 <<< nBlocksC, nThreads >>>(Nx,Ny,psihat_d,kx_d,temp1_d);
	cutilCheckMsg("Kernel execution failed: [ nextstep4 ]");
	cufftSafeCall(cufftMyExecB(planz2d,(cufftCmplx *)temp1_d,(cufftReal_ *)u_d));


	//!!!!!!!!!!!!!!!!!!!!!!!!
	//!initial nonlinear term!
	//!!!!!!!!!!!!!!!!!!!!!!!!'
	nonlin1 <<< nBlocksC, nThreads >>>(Nx,kx_d,omeghat_d,temp1_d);
	cutilCheckMsg("Kernel execution failed: [ nonlin1 ]");
    cufftSafeCall(cufftMyExecB(planz2d,(cufftCmplx *)temp1_d,(cufftReal_ *)temp2_d));
	nonlin2 <<< nBlocksR, nThreads >>>(u_d,temp2_d,nl_d);
	cutilCheckMsg("Kernel execution failed: [ nonlin2 ]");
    nonlin3 <<< nBlocksC, nThreads >>>(Ny,ky_d,omeghat_d,temp1_d);
	cutilCheckMsg("Kernel execution failed: [ nonlin3 ]");
	cufftSafeCall(cufftMyExecB(planz2d,(cufftCmplx *)temp1_d,(cufftReal_ *)temp2_d));
	nonlin4 <<< nBlocksR, nThreads >>>(Nx,Ny,nl_d,v_d,temp2_d);
	cutilCheckMsg("Kernel execution failed: [ nonlin4 ]");
    cutilSafeCall(cudaMemcpy(omeg, nl_d, Nx * Ny * sizeof(real), cudaMemcpyDeviceToHost));
    
    
	cufftSafeCall(cufftMyExecF(pland2z,(cufftReal_ *)nl_d,(cufftCmplx *)nlhat_d));
	//!!!!!!!!!!!!!!!!!!!!!

	printf("Got initial data, starting timestepping\n");	
	gettimeofday(&start_time, NULL);
	
	for(nn=1; nn<=nplots; nn++){
		for(n=1; n<=plotgap; n++){
			chg=1.0;
			copyCmplxArray <<< nBlocksC, nThreads >>>(nloldhat_d,nlhat_d);
			cutilCheckMsg("Kernel execution failed: [ copyCmplxArray ]");
			copyCmplxArray <<< nBlocksC, nThreads >>>(omegoldhat_d,omeghat_d);
			cutilCheckMsg("Kernel execution failed: [ copyCmplxArray ]");
			while(chg>tol){
				//!!!!!!!!!!!!!!!!!!!!!!
				//!{n,k} nonlinear term!
				//!!!!!!!!!!!!!!!!!!!!!!
				nonlin1 <<< nBlocksC, nThreads >>>(Ny,kx_d,omeghat_d,temp1_d);
				cutilCheckMsg("Kernel execution failed: [ nonlin1 ]");
				cufftSafeCall(cufftMyExecB(planz2d, (cufftCmplx *)temp1_d, (cufftReal_ *)temp2_d));
				nonlin2 <<< nBlocksR, nThreads >>>(u_d, temp2_d, nl_d);
				cutilCheckMsg("Kernel execution failed: [ nonlin2 ]");
				nonlin3 <<< nBlocksC, nThreads >>>(Ny, ky_d, omeghat_d, temp1_d);
				cutilCheckMsg("Kernel execution failed: [ nonlin3 ]");
				cufftSafeCall(cufftMyExecB(planz2d, (cufftCmplx *)temp1_d, (cufftReal_ *)temp2_d));
				nonlin4 <<< nBlocksR, nThreads >>>(Nx, Ny, nl_d, v_d, temp2_d);
				cutilCheckMsg("Kernel execution failed: [ nonlin4 ]");
 				cutilSafeCall(cudaMemcpy(omeg, nl_d, Nx * Ny * sizeof(real), cudaMemcpyDeviceToHost));
    			
				cufftSafeCall(cufftMyExecF(pland2z, (cufftReal_ *)nl_d, (cufftCmplx *)nlhat_d));
				//!!!!!!!!!!!!!!!!!!!!!
				nextstep1 <<< nBlocksC, nThreads >>>(Nx,Ny,dt,Re,kx_d,ky_d,omeghat_d,omegoldhat_d,nloldhat_d,nlhat_d);
				cutilCheckMsg("Kernel execution failed: [ nextstep1 ]");
				// Calculate streamfunction in fourier space, psihat
		 		nextstep2 <<< nBlocksC, nThreads >>> (Nx, Ny, kx_d, ky_d, omeghat_d, psihat_d);
				cutilCheckMsg("Kernel execution failed: [ nextstep2 ]");
 				// Calculate y velocity
				nextstep3 <<< nBlocksC, nThreads >>>(Nx,Ny,psihat_d,kx_d,temp1_d);
				cutilCheckMsg("Kernel execution failed: [ nextstep3 ]");
				cufftSafeCall(cufftMyExecB(planz2d,(cufftCmplx *)temp1_d,(cufftReal_ *)v_d));
				// Calculate x velocity
				nextstep4 <<< nBlocksC, nThreads >>>(Nx,Ny,psihat_d,kx_d,temp1_d);
				cutilCheckMsg("Kernel execution failed: [ nextstep4 ]");
				cufftSafeCall(cufftMyExecB(planz2d,(cufftCmplx *)temp1_d,(cufftReal_ *)u_d));

				cufftSafeCall(cufftMyExecB(planz2d, (cufftCmplx *)omeghat_d, (cufftReal_ *)omeg_d));
				checkConvergence1 <<< nBlocksR, nThreads >>>(Nx, Ny, omeg_d, omegcheck_d);
				cutilCheckMsg("Kernel execution failed: [ checkConvergence1 ]");
				cufftSafeCall(cufftMyExecF(pland2z, (cufftReal_ *)omegcheck_d, (cufftCmplx *)temp1_d));
				cutilSafeCall(cudaMemcpy(omeg, temp1_d, sizeof(real), cudaMemcpyDeviceToHost));
				chg=omeg[0];				
				copyRealArray <<< nBlocksR, nThreads >>>(omegcheck_d, omeg_d);
				cutilCheckMsg("Kernel execution failed: [ copyRealArray ]");	
			}
		}	
	}
	
	gettimeofday(&end_time, NULL);
	start_time.tv_sec = end_time.tv_sec - start_time.tv_sec;
	start_time.tv_usec = end_time.tv_usec - start_time.tv_usec;
	printf ("Timstepping took %lf seconds...\n", (real)(start_time.tv_sec) + (real)(start_time.tv_usec) / (real)1000000.);
	
	cutilSafeCall(cudaMemcpy(omeg, omeg_d, Nx * Ny * sizeof(real), cudaMemcpyDeviceToHost));

	// get exact omega
	for(j=0; j<Ny; j++) {
		 for(i=0;i<Nx;i++){
			omegexact[j*Ny+i]=4.0*pi*sin(2.0*pi*x[i])*sin(2.0*pi*y[j])*exp(-8.0*pi*pi*(real)nplots*(real)plotgap*dt/Re);
		}
	}
	
	max=0;
	for(i=0; i<Nx*Ny;i++) {
		omeg[i]=omeg[i]/(real)(Nx*Ny);
		chg=abs(omeg[i]-omegexact[i]);
		if(chg>=max) 
			max=chg;
	} 
	printf("Maximum error %lf ...\n", max);
	
	// turn of saving data for benchmarking
	//savedata (Nx,Ny,0,omeg);
	printf("Saved to disk\n");
	cutilSafeCall (cudaFree ((void*)kx_d));
	cutilSafeCall (cudaFree ((void*)ky_d));
	cutilSafeCall (cudaFree ((void*)x_d));
	cutilSafeCall (cudaFree ((void*)y_d));
	cutilSafeCall (cudaFree ((void*)u_d));
	cutilSafeCall (cudaFree ((void*)v_d));
	cutilSafeCall (cudaFree ((void*)temp1_d));
	cutilSafeCall (cudaFree ((void*)temp2_d));
	cutilSafeCall (cudaFree ((void*)omeg_d));
	cutilSafeCall (cudaFree ((void*)nl_d));
	cutilSafeCall (cudaFree ((void*)omegcheck_d));
	cutilSafeCall (cudaFree ((void*)omegoldhat_d));
	cutilSafeCall (cudaFree ((void*)nloldhat_d));
	cutilSafeCall (cudaFree ((void*)nlhat_d));
	cutilSafeCall (cudaFree ((void*)psihat_d));
	printf("Deallocated GPU arrays \n");

	cufftSafeCall (cufftDestroy (pland2z));
	cufftSafeCall (cufftDestroy (planz2d));
	printf("Destroyed CUFFT plans \n");
	free (kx);
	free (ky);
	free (x);
	free (y);
	free (omeg);
	free (omegexact);
	printf ("Deallocated CPU arrays\n");
	printf("End Program \n");
	return 0;
}
