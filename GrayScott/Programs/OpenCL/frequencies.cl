#if defined(cl_amd_fp64) || defined(cl_khr_fp64)
    #if defined(cl_amd_fp64)
    #pragma OPENCL EXTENSION cl_amd_fp64 : enable
    #elif defined(cl_khr_fp64)
    #pragma OPENCL EXTENSION cl_khr_fp64 : enable
    #endif
    // function declarations/definitions using double precision doubleing-point arithmetic
#endif
__kernel void frequencies ( __global double* kx, const double Lx, const int Nx)
{
   	const int ind = get_global_id(0);
	if ( ind < Nx/2) kx[ind] = -1.0*((double)((ind))/Lx)*((double)( ind)/Lx);
	if ( ind ==Nx/2) kx[ ind]=0.0;
	if ( ind > Nx/2) kx[ ind]=-1.0*(double)(Nx-ind)/Lx*(double)(Nx-ind)/Lx;
}
