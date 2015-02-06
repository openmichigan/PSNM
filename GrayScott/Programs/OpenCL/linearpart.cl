#if defined(cl_amd_fp64) || defined(cl_khr_fp64)
    #if defined(cl_amd_fp64)
    #pragma OPENCL EXTENSION cl_amd_fp64 : enable
    #elif defined(cl_khr_fp64)
    #pragma OPENCL EXTENSION cl_khr_fp64 : enable
    #endif
    // function declarations/definitions using double precision doubleing-point arithmetic
#endif
__kernel void linearpart ( __global double2* uhat, __global double2* vhat, __global const double* Kx, __global const double* Ky, const double Du, const double Dv, const double A, const double B,  const double dt, const double c, const int Nx, const int Ny)
{
   const int ind = get_global_id(0);

int i,j,k;

j=floor((double)(ind)/(double)Nx);
i=ind-j*Nx;
double uexp;
double vexp;

uexp=dt*c*(-1.0*A+Du*(Kx[i]+Ky[j]));
vexp=dt*c*(-1.0*B+Dv*(Kx[i]+Ky[j]));

uhat[ind].x=exp(uexp)*uhat[ind].x;
uhat[ind].y=exp(uexp)*uhat[ind].y;

vhat[ind].x=exp(vexp)*vhat[ind].x;
vhat[ind].y=exp(vexp)*vhat[ind].y;
}
