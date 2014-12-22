#if defined(cl_amd_fp64) || defined(cl_khr_fp64)
    #if defined(cl_amd_fp64)
    #pragma OPENCL EXTENSION cl_amd_fp64 : enable
    #elif defined(cl_khr_fp64)
    #pragma OPENCL EXTENSION cl_khr_fp64 : enable
    #endif
    // function declarations/definitions using double precision doubleing-point arithmetic
#endif
__kernel void grid ( __global double* x, const double Lx, const int Nx)
{
  const int ind = get_global_id(0);
	x[ind]=(-1.0 + ((double) 2.0*ind/(double)Nx))*M_PI*Lx;   	
	
}
