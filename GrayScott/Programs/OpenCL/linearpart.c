#if defined(cl_amd_fp64) || defined(cl_khr_fp64)
    #if defined(cl_amd_fp64)
    #pragma OPENCL EXTENSION cl_amd_fp64 : enable
    #elif defined(cl_khr_fp64)
    #pragma OPENCL EXTENSION cl_khr_fp64 : enable
    #endif
    // function declarations/definitions using double precision doubleing-point arithmetic
#endif
__kernel void linearpart ( __global double* uhatre,__global double* uhatim,__global double* vhatre,__global double* vhatim,__global const double* Kx,__global const double* Ky,__global const double* Kz, const double dt, const double Du,const double Dv, const double A,const double B,const double bhighre,const double bhighim , const int Nx, const int Ny, const int Nz)
{
   const int ind = get_global_id(0);

int i,j,k;
k=floor((double)ind/(double)(Ny*Nx));
j=floor((double)(ind-k*(Ny*Nx))/(double)Nx);
i=ind-k*(Ny*Nx)-j*Nx;
double uexp[2];
uexp[0]=dt*bhighre*(-1.0*A+Du*(Kx[i]+Ky[j]+Kz[k]));
uexp[1]=dt*bhighim*(-1.0*A+Du*(Kx[i]+Ky[j]+Kz[k]));
double vexp[2];
vexp[0]=dt*bhighre*(-1.0*B+Dv*(Kx[i]+Ky[j]+Kz[k]));
vexp[1]=dt*bhighim*(-1.0*B+Dv*(Kx[i]+Ky[j]+Kz[k]));
if(ind==0){
double N=(double)Nx*Ny*Nz;

uhatre[ind]=exp(uexp[0])*(((uhatre[ind]-N)*cos(uexp[1]))-(uhatim[ind]*sin(uexp[1])))+N;
uhatim[ind]=exp(uexp[0])*(((uhatre[ind]-N)*sin(uexp[1]))+(uhatim[ind]*cos(uexp[1])));
}

else{
uhatre[ind]=exp(uexp[0])*((uhatre[ind]*cos(uexp[1]))-(uhatim[ind]*sin(uexp[1])));
uhatim[ind]=exp(uexp[0])*((uhatre[ind]*sin(uexp[1]))+(uhatim[ind]*cos(uexp[1])));
}
vhatre[ind]=exp(vexp[0])*((vhatre[ind]*cos(vexp[1]))-(vhatim[ind]*sin(vexp[1])));
vhatim[ind]=exp(vexp[0])*((vhatre[ind]*sin(vexp[1]))+(vhatim[ind]*cos(vexp[1])));
}
