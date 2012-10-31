/*--------------------------------------------------------------------------
*
*
* PURPOSE
*
* A CUDA program to solve the nonlinear Schrodinger equation in 2 dimensions
* i*u_t+Es*|u|^2u+u_{xx}+u_{yy}=0
* using a second order time spectral splitting scheme. The program is based
* on a similar program to solve the Sine Gordon equation by P. Rigge
*
* The boundary conditions are u(x=0,y)=u(x=2*L*\pi,y) 
* and u(x,y=0)=u(x,y=2*L*\pi)
* The initial condition is u=exp(-x^2-y^2)
*
* AUTHORS
*
* B. Cloutier, B.K. Muite, P. Rigge
* 4 June 2012
*
* .. Scalars ..
*
*  Nx				= number of modes in x - power of 2 for FFT
*  Ny				= number of modes in y - power of 2 for FFT
*  dt				= timestep
*  Nt				= number of timesteps to take
*  plotgap			= number of timesteps between plots
*  Lx				= width of box in x direction
*  Ly				= width of box in y direction
*  ES				= +1 for focusing and -1 for defocusing
*  i				= loop counter in x direction
*  j				= loop counter in y direction
*  n				= loop counter for timesteps direction	
*  nThreads			= Number of threads for GPU to use
*  nBlocks			= number of blocks for GPU to use
*  plan				= fft plan
*  dt				= timestep
*  InMass			= initial mass
*  FiMass			= final mass
*  InEner			= initial energy
*  FiEner			= final energy
*  scalemodes		= scaled an array after performing inverse FFT
*  plan				= plan for fft
*  xsize    		= size of real arrays in x direction
*  ysize    		= size of real arrays in y direction
*  gridsize 		= size of array for complex data
*  complxsize 		= size of a complex data point 
*  start_time		= start time of computation
*  end_time			= end time of evaluation
*
* .. Arrays on CPU ..
*
*  u				= approximate solution
*  kx				= wave numbers in x direction
*  ky				= wave numbers in y direction
*
* .. Arrays on GPU ..
*
*  u_d				= approximate solution on device
*  v_d				= Fourier transform of approximate solution on device
*  temp1_d			= temporary array used to find mass and energy
*  temp2_d			= temporary array used to find mass and energy
*
* REFERENCES
*
* ACKNOWLEDGEMENTS
*
* ACCURACY
*		
* ERROR INDICATORS AND WARNINGS
*
* FURTHER COMMENTS
*
* Check that the initial condition is consistent with the 
* boundary conditions for the domain specified
*
* For consistency with Fortran programs, real is the same as double
* and complx is double2. The relevant complex arithmetic operations
* have been defined appropriately.
*
*--------------------------------------------------------------------
*
* External routines required
* 
* External libraries required
* cufft	 -- Cuda FFT library
*
*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <sys/time.h>

#include <cufft.h>
#include <cutil_inline.h>

typedef double2 cmplx;
typedef double real;


static __device__ __host__ inline cmplx cmul (cmplx a, cmplx b) {
  cmplx c;
  c.x = a.x * b.x - a.y * b.y;
  c.y = a.x * b.y + a.y * b.x;
  return c;
}

static __device__ inline cmplx cscale (cmplx a, real b) {
  cmplx c;
  c.x = a.x * b;
  c.y = a.y * b;
  return c;
}

static __device__ inline cmplx cexp (cmplx a) {
  cmplx c;
  c.x = exp(a.x) * cos(a.y);
  c.y = exp(a.x) * sin(a.y);
  return c;
}

static __device__ inline real abssq (cmplx a) {
  real c;
  c = (a.x)*(a.x)+(a.y)*(a.y);
  return c;
}

static __global__ void potentialcal(cmplx *v_d, cmplx *u_d, real scalemodes, real Es) 
{
  int idx = blockIdx.x * blockDim.x + threadIdx.x;
  real temp;
  temp = abssq(u_d[idx]);
  v_d[idx].x = Es*temp*scalemodes*scalemodes;
  v_d[idx].y = 0;
}

static __global__ void uxencalc(cmplx *v_d, real *kx_d, cmplx *temp1_d,
                                real scalemodes, int Ny) 
{
  int idx = blockIdx.x * blockDim.x + threadIdx.x;
  int i;
  i = idx / Ny;
  cmplx wavenum;
  wavenum.x = 0;
  wavenum.y = kx_d[i];
  cmplx temp;
  temp = cmul(wavenum,v_d[idx]);
  temp1_d[idx]=cscale(temp,scalemodes);
}

static __global__ void uyencalc(cmplx *v_d, real *ky_d, cmplx *temp1_d,
                                real scalemodes, int Ny) 
{
  int idx = blockIdx.x * blockDim.x + threadIdx.x;
  int j;
  j = idx % Ny; 
  cmplx wavenum;
  wavenum.x = 0;
  wavenum.y = ky_d[j];
  cmplx temp;
  temp = cmul(wavenum,v_d[idx]);
  temp1_d[idx]=cscale(temp,scalemodes);
}

static __global__ void potencalc(cmplx *u_d, cmplx *temp1_d, real Es ) 
{
  int idx = blockIdx.x * blockDim.x + threadIdx.x;
  real temp;
  temp = -Es*(real)0.25*abssq(u_d[idx])*abssq(u_d[idx]);
  temp1_d[idx].x = temp;
  temp1_d[idx].y = 0;
}

static __global__ void abscalc(cmplx *u_d, cmplx *temp1_d ) 
{
  int idx = blockIdx.x * blockDim.x + threadIdx.x;
  real temp;
  temp = abssq(u_d[idx]);
  temp1_d[idx].x = temp;
  temp1_d[idx].y = 0;
}
  
static __global__ void realstep(cmplx *v_d, cmplx *u_d, real scalemodes, real dt) 
{
  int idx = blockIdx.x * blockDim.x + threadIdx.x;
  cmplx imag;
  imag.x = 0;
  imag.y = 1;
  cmplx temp1;
  temp1 = cscale(v_d[idx],dt);
  cmplx temp2;
  temp2 = cmul(imag, temp1);
  cmplx temp3;
  temp3 = cexp(temp2);
  cmplx temp4;
  temp4 = cmul(temp3,u_d[idx]);
  u_d[idx] = cscale(temp4,scalemodes);
}

static __global__ void fourierstep(real *kx_d, real *ky_d,
                   cmplx *v_d, int Nx, int Ny, real dt) 
{
  int idx = blockIdx.x * blockDim.x + threadIdx.x;
  cmplx imag;
  imag.x = 0;
  imag.y = 1;
  int i, j;
  i = idx / Ny;
  j = idx % Ny; 
  real wavenum;
  wavenum =  dt*(-kx_d[i]*kx_d[i] + -ky_d[j]*ky_d[j]);
  cmplx intfactor;
  intfactor = cexp(cscale(imag,wavenum));
  v_d[j*Nx+i]=cmul(intfactor,v_d[j*Nx+i]);
 }

static __global__ void fourierstephalf(real *kx_d, real *ky_d,
                   			cmplx *v_d, int Nx, int Ny, real dt) 
{
  int idx = blockIdx.x * blockDim.x + threadIdx.x;
  cmplx imag;
  imag.x = 0;
  imag.y = 1;
  int i, j;
  i = idx / Ny;
  j = idx % Ny; 
  real wavenum;
  wavenum =  (real)0.5*dt*(-kx_d[i]*kx_d[i] + -ky_d[j]*ky_d[j]);
  cmplx intfactor;
  intfactor = cexp(cscale(imag,wavenum));
  v_d[j*Nx+i]=cmul(intfactor,v_d[j*Nx+i]);
 }

static __global__ void scalefinal(cmplx *u_d, real scalemodes) 
{
  int idx = blockIdx.x * blockDim.x + threadIdx.x;
  u_d[idx] = cscale(u_d[idx],scalemodes);
}

void initialdata (int Nx, int Ny, real *x, real *y, cmplx *u) {
  int i, j;
  for (j=0; j<Ny; j++) {
    for (i=0; i<Nx; i++) {
      u[j*Nx+i].x = 0.5*exp(-1.0*(x[i]*x[i] + y[j]*y[j]));
    }
  }
}

void savedata (int Nx, int Ny, int nplot, cmplx *u) {
  FILE *f = NULL;
  char nameconfig1[128];
  nameconfig1[0]='\0';
  sprintf (nameconfig1, "data/u%.10d.datbin", nplot);
  f = fopen (nameconfig1, "wb");
  real *ureal;
  ureal = (real*)malloc (Nx * Ny * sizeof(real));
  int i;
  for (i=0; i<Nx*Ny; i++) {
    ureal[i] = u[i].x;
  }
  fwrite (ureal, sizeof(real), Nx*Ny, f);
  fclose (f);
}

int main (int argc, char** argv) {
  cufftHandle plan;
  int Nx=128;
  int Ny=128;
  int Nt=20;
  int plotgap=1;
  int Lx=5.0;
  int Ly=5.0;
  real Es=1.0;
  real dt=1.e-5;
  real *kx, *ky, *kx_d, *ky_d, *x, *y;
  cmplx *u, *u_d, *v_d, *temp1_d, *temp2_d;
  size_t xsize    = Nx * sizeof(real);
  size_t ysize    = Ny * sizeof(real);
  size_t gridsize = Nx * Ny * sizeof(cmplx);
  size_t complxsize = sizeof(cmplx);
  struct timeval start_time, end_time;
  int nThreads, nBlocks;
  int i, j, n;
  real scalemodes;
  cmplx InMass, FiMass, InEner, FiEner;

  printf("Program starting\n");
  printf("Grid: %d X %d\n",Nx,Ny);
  printf("dt: %lf\n",dt);
  
  scalemodes=(real)1.0 / ( ( (real)Nx) * ( (real)Ny) );
  nThreads = 256;
  nBlocks  = Nx * Ny / nThreads;

  /* Allocate */
  kx = (real*) malloc (xsize);
  ky = (real*) malloc (ysize);
  x  = (real*) malloc (xsize);
  y  = (real*) malloc (ysize);
  u   = (cmplx*) malloc (gridsize);

  /* Plan FFTs */
  cufftSafeCall (
      cufftPlan2d (&plan, Nx, Ny, CUFFT_Z2Z));
  printf ("Set up FFTs...\n");

  cutilSafeCall (
      cudaMalloc ((void**)&kx_d, xsize));
  cutilSafeCall (
      cudaMalloc ((void**)&ky_d, ysize));
  cutilSafeCall (
      cudaMalloc ((void**)&u_d,  gridsize));
  cutilSafeCall (
      cudaMalloc ((void**)&v_d,  gridsize));
  cutilSafeCall (
      cudaMalloc ((void**)&temp1_d,  gridsize));
  cutilSafeCall (
      cudaMalloc ((void**)&temp2_d,  gridsize));
  printf ("Allocated GPU arrays...\n");

  /* Initialize arrays */
  for (i=0; i<Nx/2; i++) {
    kx[i] = i/(real)Lx;
  }
  kx[Nx/2]=0.;
  for (i=0; i<Nx/2-1; i++) {
    kx[Nx/2+1+i] = -kx[Nx/2-i-1];
  }
  for (i=0; i<Ny/2; i++) {
    ky[i] = i/(real)Ly;
  }
  ky[Ny/2]=0.;
  for (i=0; i<Ny/2-1; i++) {
    ky[Ny/2+1+i] = -ky[Ny/2-i-1];
  }
  for (i=0; i<Nx; i++) {
    x[i] = (-1. + (2. * i)/(real)Nx) * M_PI * Lx;
  }
  for (j=0; j<Ny; j++) {
    y[j] = (-1. + (2. * j)/(real)Ny) * M_PI * Ly;
  }
  printf ("Initialized arrays...\n");

  /* Initial data */
  initialdata (Nx, Ny, x, y, u);
  printf ("Got initial data...\n");
  
  savedata(Nx, Ny, 0, u);
  printf ("Saved initial data...\n");

  cutilSafeCall (
      cudaMemcpy (kx_d, kx, xsize, cudaMemcpyHostToDevice));
  cutilSafeCall (
      cudaMemcpy (ky_d, ky, ysize, cudaMemcpyHostToDevice));
  cutilSafeCall (
    cudaMemcpy (u_d,  u, gridsize, cudaMemcpyHostToDevice));
  printf ("Copied initial data to device...\n");
  cufftSafeCall (
      cufftExecZ2Z (plan,
                    (cufftDoubleComplex *)u_d,
                    (cufftDoubleComplex *)v_d,
                    CUFFT_FORWARD));

  printf ("Arrays initialized...\n");
  /* Calculate initial mass */
  abscalc <<< nBlocks, nThreads >>> (u_d, temp1_d);
  cutilCheckMsg("Kernel execution failed: [ abscalc ]");
    cufftSafeCall (
        cufftExecZ2Z (plan,
                      (cufftDoubleComplex *)temp1_d,
                      (cufftDoubleComplex *)temp2_d,
                      CUFFT_FORWARD));
  cutilSafeCall (
      cudaMemcpy (u, temp2_d, complxsize, cudaMemcpyDeviceToHost));
   InMass=u[0];
  /* Calculate initial energy */
  potencalc <<< nBlocks, nThreads >>> (u_d, temp1_d, Es);    
  cutilCheckMsg("Kernel execution failed: [ potencalc ]");
    cufftSafeCall (
        cufftExecZ2Z (plan,
                      (cufftDoubleComplex *)temp1_d,
                      (cufftDoubleComplex *)temp2_d,
                      CUFFT_FORWARD));
  cutilSafeCall (
      cudaMemcpy (u, temp2_d, complxsize, cudaMemcpyDeviceToHost));
  InEner=u[0];
  uxencalc <<< nBlocks, nThreads >>> (v_d, kx_d, temp1_d, scalemodes, Ny);    
  cutilCheckMsg("Kernel execution failed: [ uxencalc ]");
    cufftSafeCall (
        cufftExecZ2Z (plan,
                      (cufftDoubleComplex *)temp1_d,
                      (cufftDoubleComplex *)temp2_d,
                      CUFFT_INVERSE));
  abscalc <<< nBlocks, nThreads >>> (temp2_d, temp1_d);
  cutilCheckMsg("Kernel execution failed: [ abscalc ]");
    cufftSafeCall (
        cufftExecZ2Z (plan,
                      (cufftDoubleComplex *)temp2_d,
                      (cufftDoubleComplex *)temp1_d,
                      CUFFT_FORWARD));
  cutilSafeCall (
      cudaMemcpy (u, temp1_d, complxsize, cudaMemcpyDeviceToHost));
   InEner.x=InEner.x+(real)0.5*u[0].x;
   InEner.y=InEner.y+(real)0.5*u[0].y;
  uyencalc <<< nBlocks, nThreads >>> (v_d, ky_d, temp1_d, scalemodes, Ny);    
  cutilCheckMsg("Kernel execution failed: [ uyencalc ]");
    cufftSafeCall (
        cufftExecZ2Z (plan,
                      (cufftDoubleComplex *)temp1_d,
                      (cufftDoubleComplex *)temp2_d,
                      CUFFT_INVERSE));
  abscalc <<< nBlocks, nThreads >>> (temp2_d, temp1_d);
  cutilCheckMsg("Kernel execution failed: [ abscalc ]");
  cufftSafeCall (
        cufftExecZ2Z (plan,
                      (cufftDoubleComplex *)temp2_d,
                      (cufftDoubleComplex *)temp1_d,
                      CUFFT_FORWARD));
  cutilSafeCall (
      cudaMemcpy (u, temp1_d, complxsize, cudaMemcpyDeviceToHost));
   InEner.x=InEner.x+(real)0.5*u[0].x;
   InEner.y=InEner.y+(real)0.5*u[0].y;
  
                
  gettimeofday(&start_time, NULL);
  /* Do first half time step */	
  fourierstephalf <<< nBlocks, nThreads >>> (kx_d, ky_d, v_d, Nx, Ny, dt);
  cutilCheckMsg("Kernel execution failed: [ fourierstephalf ]");
   

  for (n=0; n<Nt; n++) {
    cufftSafeCall (
        cufftExecZ2Z (plan,
                      (cufftDoubleComplex *)v_d,
                      (cufftDoubleComplex *)u_d,
                      CUFFT_INVERSE));
    potentialcal <<< nBlocks, nThreads >>>(v_d, u_d, scalemodes, Es);
    cutilCheckMsg("Kernel execution failed: [potentialcal]");
    realstep <<< nBlocks, nThreads >>>(v_d, u_d, scalemodes, dt);
    cutilCheckMsg("Kernel execution failed: [realstep]");
    cufftSafeCall (
        cufftExecZ2Z (plan,
                      (cufftDoubleComplex *)u_d,
                      (cufftDoubleComplex *)v_d,
                      CUFFT_FORWARD));
    fourierstep <<< nBlocks, nThreads >>> (kx_d, ky_d, v_d, Nx, Ny, dt);
    cutilCheckMsg("Kernel execution failed: [ fourierstep ]");
  }
  
  /* transform back final data and do another half step */
  cufftSafeCall (
        cufftExecZ2Z (plan,
                      (cufftDoubleComplex *)v_d,
                      (cufftDoubleComplex *)u_d,
                      CUFFT_INVERSE));
  potentialcal <<< nBlocks, nThreads >>>(v_d, u_d, scalemodes, Es);
  cutilCheckMsg("Kernel execution failed: [potentialcal]");
  realstep <<< nBlocks, nThreads >>>(v_d, u_d, scalemodes, dt);
  cutilCheckMsg("Kernel execution failed: [realstep]");
  cufftSafeCall (
        cufftExecZ2Z (plan,
                      (cufftDoubleComplex *)u_d,
                      (cufftDoubleComplex *)v_d,
                      CUFFT_FORWARD));
  
  fourierstephalf <<< nBlocks, nThreads >>> (kx_d, ky_d, v_d, Nx, Ny, dt);
  cutilCheckMsg("Kernel execution failed: [ fourierstephalf ]");
  
  cufftSafeCall (
        cufftExecZ2Z (plan,
                      (cufftDoubleComplex *)v_d,
                      (cufftDoubleComplex *)u_d,
                      CUFFT_INVERSE));
  scalefinal <<< nBlocks, nThreads >>>(u_d, scalemodes);
  
  gettimeofday(&end_time, NULL);
  start_time.tv_sec = end_time.tv_sec - start_time.tv_sec;
  start_time.tv_usec = end_time.tv_usec - start_time.tv_usec;

  printf ("Computation complete...\n");
  printf ("Computation took %lf seconds...\n",
  (real)(start_time.tv_sec) + (real)(start_time.tv_usec) / (real)1000000);
  cutilSafeCall (
          cudaMemcpy ( u, u_d, gridsize, cudaMemcpyDeviceToHost));
  savedata(Nx,Ny,1+n/plotgap, u); 
  /* Calculate final mass */
  abscalc <<< nBlocks, nThreads >>> (u_d, temp1_d);
  cutilCheckMsg("Kernel execution failed: [ abscalc ]");
    cufftSafeCall (
        cufftExecZ2Z (plan,
                      (cufftDoubleComplex *)temp1_d,
                      (cufftDoubleComplex *)temp2_d,
                      CUFFT_FORWARD));
  cutilSafeCall (
      cudaMemcpy (u, temp2_d, complxsize, cudaMemcpyDeviceToHost));
 FiMass=u[0];
   /* Calculate Final energy */
  potencalc <<< nBlocks, nThreads >>> (u_d, temp1_d, Es);    
  cutilCheckMsg("Kernel execution failed: [ potencalc ]");
    cufftSafeCall (
        cufftExecZ2Z (plan,
                      (cufftDoubleComplex *)temp1_d,
                      (cufftDoubleComplex *)temp2_d,
                      CUFFT_FORWARD));
  cutilSafeCall (
      cudaMemcpy (u, temp2_d, complxsize, cudaMemcpyDeviceToHost));
  FiEner=u[0];
  uxencalc <<< nBlocks, nThreads >>> (v_d, kx_d, temp1_d, scalemodes, Ny);    
  cutilCheckMsg("Kernel execution failed: [ uxencalc ]");
    cufftSafeCall (
        cufftExecZ2Z (plan,
                      (cufftDoubleComplex *)temp1_d,
                      (cufftDoubleComplex *)temp2_d,
                      CUFFT_INVERSE));
  abscalc <<< nBlocks, nThreads >>> (temp2_d, temp1_d);
  cutilCheckMsg("Kernel execution failed: [ abscalc ]");
    cufftSafeCall (
        cufftExecZ2Z (plan,
                      (cufftDoubleComplex *)temp2_d,
                      (cufftDoubleComplex *)temp1_d,
                      CUFFT_FORWARD));
  cutilSafeCall (
      cudaMemcpy (u, temp1_d, complxsize, cudaMemcpyDeviceToHost));
   FiEner.x=FiEner.x+(real)0.5*u[0].x;
   FiEner.y=FiEner.y+(real)0.5*u[0].y;
  uyencalc <<< nBlocks, nThreads >>> (v_d, ky_d, temp1_d, scalemodes, Ny);    
  cutilCheckMsg("Kernel execution failed: [ uyencalc ]");
    cufftSafeCall (
        cufftExecZ2Z (plan,
                      (cufftDoubleComplex *)temp1_d,
                      (cufftDoubleComplex *)temp2_d,
                      CUFFT_INVERSE));
  abscalc <<< nBlocks, nThreads >>> (temp2_d, temp1_d);
  cutilCheckMsg("Kernel execution failed: [ abscalc ]");
    cufftSafeCall (
        cufftExecZ2Z (plan,
                      (cufftDoubleComplex *)temp2_d,
                      (cufftDoubleComplex *)temp1_d,
                      CUFFT_FORWARD));
  cutilSafeCall (
      cudaMemcpy (u, temp1_d, complxsize, cudaMemcpyDeviceToHost));
   FiEner.x=FiEner.x+(real)0.5*u[0].x;
   FiEner.y=FiEner.y+(real)0.5*u[0].y;

 printf ("Initial mass %lf ...\n", InMass.x);
 printf ("Final mass %lf ...\n", FiMass.x);
 printf ("Initial energy %lf ...\n", InEner.x);
 printf ("Final energy %lf ...\n", FiEner.x);
 
 cutilSafeCall (
      cudaFree ((void*)kx_d));
  cutilSafeCall (
      cudaFree ((void*)ky_d));
  cutilSafeCall (
      cudaFree ((void*)u_d));
  cutilSafeCall (
      cudaFree ((void*)v_d));
  cutilSafeCall (
      cudaFree ((void*)temp1_d));
  cutilSafeCall (
      cudaFree ((void*)temp2_d));
  cufftSafeCall (
      cufftDestroy (plan));
  free (kx);
  free (ky);
  free (x);
  free (y);
  free (u);
  printf ("Cleaned up...\n");
  return 0;
}
