#if defined(cl_amd_fp64) || defined(cl_khr_fp64)
    #if defined(cl_amd_fp64)
    #pragma OPENCL EXTENSION cl_amd_fp64 : enable
    #elif defined(cl_khr_fp64)
    #pragma OPENCL EXTENSION cl_khr_fp64 : enable
    #endif
    // function declarations/definitions using double precision doubleing-point arithmetic
#endif
__kernel void nonlinearpart_a ( __global double2* u,__global double2* v, const double A, const double dt, const double a){
    const int ind = get_global_id(0);
	//u[ind].x=A/(v[ind].x*v[ind].x)+(u[ind].x-A/(v[ind].x*v[ind].x))*exp(dt*a*(-1.0)*v[ind].x*v[ind].x);
	u[ind].x=u[ind].x*exp(dt*a*(-1.0)*v[ind].x*v[ind].x);	
	u[ind].y=0.0;
}
