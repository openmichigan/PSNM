#if defined(cl_amd_fp64) || defined(cl_khr_fp64)
    #if defined(cl_amd_fp64)
    #pragma OPENCL EXTENSION cl_amd_fp64 : enable
    #elif defined(cl_khr_fp64)
    #pragma OPENCL EXTENSION cl_khr_fp64 : enable
    #endif
    // function declarations/definitions using double precision doubleing-point arithmetic
#endif
__kernel void initialdata ( __global double* ure, __global double* vre, __global double* uim, __global double* vim, __global const double* x ,__global const double* y ,__global const double* z , const int Nx, const int Ny, const int Nz)
{
   const int ind = get_global_id(0);

int i,j,k;
k=floor((double)ind/(double)(Ny*Nx));
j=floor((double)(ind-k*(Ny*Nx))/(double)Nx);
i=ind-k*(Ny*Nx)-j*Nx;
ure[ind]=0.5+exp(-1.0*(x[i]*x[i]+y[j]*y[j]+z[k]*z[k] )-1.0);// 	
vre[ind]=0.1+exp(-1.0*(x[i]*x[i]+y[j]*y[j]+z[k]*z[k] )-1.0);
uim[ind]=0.0;
vim[ind]=0.0;
}
