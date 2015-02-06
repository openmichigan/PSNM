#if defined(cl_amd_fp64) || defined(cl_khr_fp64)
    #if defined(cl_amd_fp64)
    #pragma OPENCL EXTENSION cl_amd_fp64 : enable
    #elif defined(cl_khr_fp64)
    #pragma OPENCL EXTENSION cl_khr_fp64 : enable
    #endif
    // function declarations/definitions using double precision doubleing-point arithmetic
#endif
__kernel void nonlinearpart_b ( __global double2* u,__global double2* v, const double A, const double dt, const double b){
    const int ind = get_global_id(0);
	//v[ind].x=-v[ind].x/(u[ind].x*dt*b*v[ind].x-1.0);

	v[ind].x=1.0/(-0.5*A*(dt*dt*b*b)-u[ind].x*(dt*b)+1/v[ind].x);
	v[ind].y=0.0;
	u[ind].x=A*dt*b+u[ind].x;
}
