/*--------------------------------------------------------------------
*
*
* PURPOSE
*
* This program solves nonlinear sine-Gordon equation in 2 dimensions
* u_{tt}-u_{xx}-u_{yy}=-sin(u)
* using a second order implicit-explicit time stepping scheme.
*
* The boundary conditions are u(x=0,y)=u(2*Lx*\pi,y),
*       u(x,y=0)=u(x,y=2*Ly*\pi)
* 
* AUTHORS
*
* B. Cloutier, B.K. Muite, P. Rigge
* 4 June 2012
*
* .. Scalars ..
*  Nx                           = number of modes in x - power of 2 for FFT
*  Ny                           = number of modes in y - power of 2 for FFT
*  Nt                           = number of timesteps to take
*  plotgap                      = number of timesteps between plots
*  Lx                           = width of box in x direction
*  Ly                           = width of box in y direction
*  i                            = loop counter in x direction
*  j                            = loop counter in y direction
*  n                            = loop counter for timesteps 
*  nThreads						= Number of threads for GPU to use
*  nBlocksR						= number of blocks for GPU to use for real arrays
*  nBlocksC						= number of blocks for GPU to use for complex arrays
*  planfc                       = Forward 2d fft plan  (FFTW)
*  planbc                       = Backward 2d fft plan (FFTW)
*  planf                        = Forward 2d fft plan  (CUFFT)
*  planb                        = Backward 2d fft plan (CUFFT)
*  dt                           = timestep
*  xsize    					= size of real arrays in x direction
*  ysize    					= size of real arrays in y direction
*  gridsize 					= size of array 
*  start_time					= start time of computation
*  end_time						= end time of evaluation
*  en                           = total energy
*  es       	                = strain energy
*  ep     	                    = potential energy
*  ek	                        = kinetic energy
*
* .. Arrays on CPU ..
*
*  u                            = approximate solution
*  uold                         = approximate solution at previous timestep
*  temp1                        = extra space for energy computation
*  temp2                        = extra space for energy computation
*  kx                           = fourier frequencies in x direction (real format)
*  ky                           = fourier frequencies in y direction (real format)
*  kx_c                         = fourier frequencies in x direction (complex format)
*  ky_c                         = fourier frequencies in y direction (complex format)
*  x                            = x locations
*  y                            = y locations
*
* .. Arrays on GPU ..
*
*  u_d							= approximate solution
*  v_d							= Fourier transform of approximate solution
*  vold_d						= Fourier transform of approximate solution at previous timestep
*  nonlinhat_d					= Fourier transform of nonlinear term
*  kx_d                         = fourier frequencies in x direction (on GPU)
*  ky_d                         = fourier frequencies in y direction (on GPU)
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
* Check that the initial iterate is consistent with the
* boundary conditions for the domain specified
*
*--------------------------------------------------------------------
*
* External routines required
*       getgrid.f90     -- Get initial grid of points
*       initialdata.f90 -- Get initial data
*       enercalc.f90    -- Subroutine to calculate the energy
*       savedata.f90    -- Save initial data
*
* External libraries required
*       Cuda FFT        -- http://developer.nvidia.com/cufft
*       FFTW3           -- Fastest Fourier Transform in the West
*                       (http://www.fftw.org/)
*       OpenMP
*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>

#include <cufft.h>
#include <cutil_inline.h>
#include <fftw3.h>


typedef double2 cmplx;
typedef double real;
typedef cufftDoubleReal cufftReal_;
typedef cufftDoubleComplex cufftCmplx;
#define cufftMyExecF cufftExecD2Z
#define cufftMyExecB cufftExecZ2D
#define CUFFT_TFORM_FORWARD  CUFFT_D2Z
#define CUFFT_TFORM_BACKWARD CUFFT_Z2D


#define IDX(i,j,Nx) (Nx*(j)+(i))


extern "C" {
extern int enercalc(int*, int*, fftw_plan*, fftw_plan*,
                    double*, double*, double*, double*,
                    double*, double**, double**, cmplx**,
                    cmplx**, double**, double**);
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

static __device__ inline cmplx cscale (cmplx a, real b) {
  cmplx c;
  c.x = a.x * b;
  c.y = a.y * b;
  return c;
}

/* Compute nonlinear term */
static __global__ void nonlinterm(real *u_d) {
  int idx = blockIdx.x * blockDim.x + threadIdx.x;
  u_d[idx] = -1.0*sin(u_d[idx]);
}

/* Scale data by 1/(Nx*Ny) after ifft */
static __global__ void scaledata(int Nx, int Ny, real *u) {
  int idx = blockIdx.x * blockDim.x + threadIdx.x;
  real s = (real)1. / ((real)Nx * (real)Ny);
  u[idx] = s * u[idx];
}

/* Compute next timestep in Fourier space */
static __global__ void nextstep(int Nx, int Ny, real dt,
								real *kx_d, real *ky_d,
                                cmplx *v_d, cmplx *vold_d,
                                cmplx *nonlinhat_d) {
  int idx = blockIdx.x * blockDim.x + threadIdx.x;
  int i, j;
  i = idx / (Ny/2+1);
  j = idx % (Ny/2+1);
  /* Minus kx and minus ky because kx, ky multiplied by i=sqrt(-1) */
  real coef1 = (real)0.25 * (-kx_d[i]*kx_d[i] - ky_d[j]*ky_d[j]);
  real coef2 = (real)1./(dt*dt);
  real coef3 = (real)1.;
  real coef4 = (real)1. / (coef2 - coef1);
  coef1 = coef1 * coef4;
  coef2 = coef2 * coef4;
  coef3 = coef3 * coef4;
  cmplx term1 = cadd (cscale (v_d[idx], (real)2.), vold_d[idx]);
  cmplx term2 = csub (cscale (v_d[idx], (real)2.), vold_d[idx]);
  cmplx term3 = nonlinhat_d[idx];
  term1 = cscale (term1, coef1);
  term2 = cscale (term2, coef2);
  term3 = cscale (term3, coef3);
  vold_d[idx] = v_d[idx];
  v_d[idx] = cadd (cadd (term1, term2), term3);
  nonlinhat_d[idx] = cscale (v_d[idx], 1./ ((real)(Nx * Ny)));
}

/* Set u and uold to initial values */
void initialdata (int Nx, int Ny, real *x, real *y, real *u, real *uold) {
  int i, j;
  for (i=0; i<Nx; i++) {
    for (j=0; j<Ny; j++) {
      u[IDX(i,j,Nx)] = (real)0.5 * exp (-x[i]*x[i]-y[j]*y[j]);
      uold[IDX(i,j,Nx)] = u[IDX(i,j,Nx)];
    }
  }
}

/* Save u in binary format  */
void savedata (int Nx, int Ny, int nplot, real *u) {
  FILE *f = NULL;
  char nameconfig1[128];
  nameconfig1[0]='\0';
  sprintf (nameconfig1, "data/u%.10d.datbin", nplot);
  f = fopen (nameconfig1, "wb");
  fwrite (u, sizeof(real), Nx*Ny, f);
  fclose (f);
}

int main (int argc, char** argv) {
  cufftHandle planf, planb;
  int Nx=1024;
  int Ny=1024;
  int Nt=500;
  int plotgap=500;
  double Lx=5.e0;
  double Ly=5.e0;
  double dt=1.e-3;
  real *kx, *ky, *kx_d, *ky_d, *x, *y, *kx_c, *ky_c;
  real *u, *u_d, *uold;
  double ek, es, ep, en;
  cmplx *vold_d, *v_d, *nonlinhat_d;;
  size_t xsize    = Nx * sizeof(real);
  size_t ysize    = Ny * sizeof(real);
  size_t gridsize = Nx * (Ny/2+1) * sizeof(cmplx);
  struct timeval start_time, end_time;
  int i, j, n;
  int nThreads, nBlocksR, nBlocksC;
  cmplx *temp1,*temp2;
  fftw_plan planfc, planbc;

  nThreads = (min(512, Nx));
  nBlocksR = (Nx * Ny / nThreads);
  nBlocksC = (Nx * (Ny/2+1) / nThreads);

  /* Print run information */
  printf ("Nx: %d\n", Nx);
  printf ("Ny: %d\n", Ny);
  printf ("Nt: %d\n", Nt);
  printf ("Lx: %lf\n", Lx);
  printf ("Ly: %lf\n", Ly);
  printf ("dt: %lf\n", dt);

  /* Allocate host arrays */
  kx  = (real*) malloc (xsize);
  ky  = (real*) malloc (ysize);
  kx_c  = (real*) malloc (xsize*2);
  ky_c  = (real*) malloc (ysize*2);
  x   = (real*) malloc (xsize);
  y   = (real*) malloc (ysize);
  u   = (real*) malloc (Nx * Ny * sizeof(real));
  uold= (real*) malloc (Nx * Ny * sizeof(real));
  temp1 = (cmplx*)malloc(Nx*Ny*sizeof(cmplx));
  temp2 = (cmplx*)malloc(Nx*Ny*sizeof(cmplx));

  /* Get cuda device */
  if (cutCheckCmdLineFlag (argc, (const char**)argv, "device")) {
    cutilDeviceInit(argc, argv);
  } else {
    cudaSetDevice (cutGetMaxGflopsDeviceId());
  }

  /* Plan cuda FFTs */
  cufftSafeCall (
      cufftPlan2d (&planf, Nx, Ny, CUFFT_TFORM_FORWARD));
  cufftSafeCall (
      cufftPlan2d (&planb, Nx, Ny, CUFFT_TFORM_BACKWARD));

  /* Plan fftw FFTs */
  planfc=fftw_plan_dft_2d(Nx,Ny,(fftw_complex*)temp1,(fftw_complex*)temp2,FFTW_FORWARD,FFTW_ESTIMATE);
  planbc=fftw_plan_dft_2d(Nx,Ny,(fftw_complex*)temp2,(fftw_complex*)temp1,FFTW_BACKWARD,FFTW_ESTIMATE);
  printf ("Set up FFTs...\n");

  /* Allocate GPU arrays */
  cutilSafeCall (
      cudaMalloc ((void**)&u_d, Nx * Ny * sizeof(cufftReal_)));
  cutilSafeCall (
      cudaMalloc ((void**)&v_d,  gridsize));
  cutilSafeCall (
      cudaMalloc ((void**)&vold_d,  gridsize));
  cutilSafeCall (
      cudaMalloc ((void**)&nonlinhat_d,  gridsize));
  cutilSafeCall (
      cudaMalloc ((void**)&kx_d, xsize));
  cutilSafeCall (
      cudaMalloc ((void**)&ky_d, ysize));
  printf ("Allocated GPU arrays...\n");

  for (i=0; i<Nx/2; i++) {
    kx[i] = i/(real)Lx;
  }
  kx[Nx/2]=0.;
  for (i=0; i<Nx/2-1; i++) {
    kx[Nx/2+1+i] = -kx[Nx/2-i-1];
  }
  for (j=0; j<Ny/2; j++) {
    ky[j] = j/(real)Ly;
  }
  ky[Ny/2]=0.;
  for (j=0; j<Ny/2-1; j++) {
    ky[Ny/2+1+j] = -1*ky[Ny/2-j-1];
  }
  for (i=0; i<Nx; i++) {
    kx_c[2*i+0]=(real)0.;
    kx_c[2*i+1]=kx[i];
  }
  for (i=0; i<Ny; i++) {
    ky_c[2*i+0]=(real)0.;
    ky_c[2*i+1]=ky[i];
  }
  for (i=0; i<Nx; i++) {
    x[i] = (-1. + (2. * i)/(real)Nx) * M_PI * Lx;
  }
  for (j=0; j<Ny; j++) {
    y[j] = (-1. + (2. * j)/(real)Ny) * M_PI * Ly;
  }

  /* Set u, uold */
  initialdata (Nx,Ny,x,y,u,uold);
  /* savedata(Nx,Ny,0, u); */ /* disabled for benchmarking */

  cutilSafeCall (
      cudaMemcpy (kx_d, kx, xsize, cudaMemcpyHostToDevice));
  cutilSafeCall (
      cudaMemcpy (ky_d, ky, ysize, cudaMemcpyHostToDevice));
  cutilSafeCall (
      cudaMemcpy (u_d, uold, Nx * Ny * sizeof(real), cudaMemcpyHostToDevice));

  cufftSafeCall (
      cufftMyExecF (planf,
                    (cufftReal_ *)u_d,
                    (cufftCmplx *)vold_d));
  cutilSafeCall (
      cudaMemcpy (u_d, u, Nx * Ny * sizeof(real), cudaMemcpyHostToDevice));

  cufftSafeCall (
      cufftMyExecF (planf,
                    (cufftReal_ *)u_d,
                    (cufftCmplx *)v_d));
  cutilSafeCall (
      cudaMemcpy (u_d, u, Nx * Ny * sizeof(real), cudaMemcpyHostToDevice));

  printf ("Arrays initialized...\n");

  enercalc(&Nx,&Ny,&planfc,&planbc,&dt,&ek,&es,&ep,&en,
            &kx_c,&ky_c,&temp1,&temp2,&u,&uold);
  printf ("Initial Energy: Tot: %lf\tKin: %lf\tStr: %lf\tPot: %lf\n", en, ek, es, ep);

  gettimeofday(&start_time, NULL);
  for (n=0; n<Nt; n++) {
    if (((n+1)%plotgap) == 0 && 0) { /* turn of plotting for tests */
      cutilSafeCall (
          cudaMemcpy (u, u_d, Nx * Ny * sizeof(real), cudaMemcpyDeviceToHost));
      /* savedata(Nx,Ny,(n+1)/plotgap, u); */ /* disabled for benchmarking */
      printf ("%lf\n", (n+1)/plotgap);
    }

    nonlinterm <<< nBlocksR, nThreads >>> (u_d);
      cutilCheckMsg("Kernel execution failed: [ nonlinterm ]");

    cufftSafeCall (
        cufftMyExecF (planf,
                      (cufftReal_ *)u_d,
                      (cufftCmplx *)nonlinhat_d));

    nextstep <<< nBlocksC, nThreads >>>(Nx, Ny, dt, kx_d, ky_d, v_d, vold_d, nonlinhat_d);
    cutilCheckMsg("Kernel execution failed: [nextstep]");
    /* scaled v_d stored in nonlinhat_d */

    cufftSafeCall (
        cufftMyExecB (planb,
                      (cufftCmplx *)nonlinhat_d,
                      (cufftReal_ *)u_d));
  }
  gettimeofday(&end_time, NULL);
  start_time.tv_sec = end_time.tv_sec - start_time.tv_sec;
  start_time.tv_usec = end_time.tv_usec - start_time.tv_usec;

  printf ("Computation complete...\n");
  printf ("Computation took %lf seconds...\n", (real)(start_time.tv_sec) + (real)(start_time.tv_usec) / (real)1000000.);

  cufftSafeCall (
      cufftMyExecB (planb,
                    (cufftCmplx *)v_d,
                    (cufftReal_ *)u_d));
  scaledata <<< nBlocksR, nThreads >>> (Nx, Ny, u_d);
  cutilCheckMsg("Kernel execution failed: [scaledata]");
  cutilSafeCall (
      cudaMemcpy (u, u_d, Nx * Ny * sizeof(real), cudaMemcpyDeviceToHost));

  cufftSafeCall (
      cufftMyExecB (planb,
                   (cufftCmplx*)vold_d,
                    (cufftReal_*)u_d));
  scaledata <<< nBlocksR, nThreads >>> (Nx,Ny,u_d);
  cutilCheckMsg("Kernel execution failed: [scaledata]");
  cutilSafeCall (
      cudaMemcpy (uold, u_d, Nx * Ny * sizeof(real), cudaMemcpyDeviceToHost));

  enercalc(&Nx,&Ny,&planfc,&planbc,&dt,&ek,&es,&ep,&en,
            &kx_c,&ky_c,&temp1,&temp2,&u,&uold);
  printf ("Final Energy: Tot: %lf\tKin: %lf\tStr: %lf\tPot: %lf\n", en, ek, es, ep);

  cutilSafeCall (
      cudaFree ((void*)kx_d));
  cutilSafeCall (
      cudaFree ((void*)ky_d));
  cutilSafeCall (
      cudaFree ((void*)u_d));
  cutilSafeCall (
      cudaFree ((void*)v_d));
  cutilSafeCall (
      cudaFree ((void*)vold_d));
  cutilSafeCall (
      cudaFree ((void*)nonlinhat_d));
  cufftSafeCall (
      cufftDestroy (planf));
  cufftSafeCall (
      cufftDestroy (planb));
  fftw_destroy_plan(planfc);
  fftw_destroy_plan(planbc);
  free (kx);
  free (ky);
  free (kx_c);
  free (ky_c);
  free (x);
  free (y);
  free (u);
  free (uold);
  free (temp1);
  free (temp2);
  printf ("Cleaned up...\n");
  return 0;
}
