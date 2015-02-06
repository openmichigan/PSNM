#if defined(cl_amd_fp64) || defined(cl_khr_fp64)
    #if defined(cl_amd_fp64)
    #pragma OPENCL EXTENSION cl_amd_fp64 : enable
    #elif defined(cl_khr_fp64)
    #pragma OPENCL EXTENSION cl_khr_fp64 : enable
    #endif
    // function declarations/definitions using double precision doubleing-point arithmetic
#endif
__kernel void initialdata ( __global double2* u, __global double2* v, const int Nx, const int Ny, const double Lx, const double Ly)
{
    const int ind = get_global_id(0);
	int i,j;
	j=floor((double)(ind)/(double)Nx);
	i=ind-j*Nx;
	double x=(-1.0 + ( 2.0*(double)i/(double)Nx))*M_PI*Lx;
	double y=(-1.0 + ( 2.0*(double)j/(double)Ny))*M_PI*Ly;
	u[ind].x=0.5+exp(-1.0*(x*x+y*y)-1.0);// 
	u[ind].y=0.0;	
	v[ind].x=0.1+exp(-1.0*(x*x+y*y)-1.0);
	v[ind].y=0.0;
}
