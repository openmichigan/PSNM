#if defined(cl_amd_fp64) || defined(cl_khr_fp64)
    #if defined(cl_amd_fp64)
    #pragma OPENCL EXTENSION cl_amd_fp64 : enable
    #elif defined(cl_khr_fp64)
    #pragma OPENCL EXTENSION cl_khr_fp64 : enable
    #endif
    // function declarations/definitions using double precision doubleing-point arithmetic
    #endif
__kernel void abserr ( __global const double* inre, __global const double* inim, __global double* outre, __global const double* outim){
   const int ind = get_global_id(0);
outre[ind]=sqrt((outre[ind]-inre[ind])*(outre[ind]-inre[ind])+(outim[ind]-inim[ind])*(outim[ind]-inim[ind]));
}
