#if defined(cl_amd_fp64) || defined(cl_khr_fp64)
    #if defined(cl_amd_fp64)
    #pragma OPENCL EXTENSION cl_amd_fp64 : enable
    #elif defined(cl_khr_fp64)
    #pragma OPENCL EXTENSION cl_khr_fp64 : enable
    #endif
    // function declarations/definitions using double precision doubleing-point arithmetic
#endif
__kernel void nonlinearpart ( __global double* ure,__global double* uim,__global double* vre,__global double* vim, const double dt, const double are,const double aim)
{
   const int ind = get_global_id(0);
	const double tol=pown(0.1,12);//0.000000000001;
double chg=1;
const double uoldre=ure[ind];
const double uoldim=uim[ind];
const double voldre=vre[ind];
const double voldim=vim[ind];
double utempre,utempim,vtempre,vtempim;
double uMre,uMim,vMre,vMim;
while(chg>tol){
utempre=ure[ind];
utempim=uim[ind];
vtempre=vre[ind];
vtempim=vim[ind];	

uMre=0.5*(ure[ind]+uoldre);
uMim=0.5*(uim[ind]+uoldim);
vMre=0.5*(vre[ind]+voldre);
vMim=0.5*(vim[ind]+voldim);								
							
ure[ind]=uoldre-dt*are*uMre*vMre*vMre+dt*are*uMre*vMim*vMim+2.0*dt*are*uMim*vMre*vMim+2.0*dt*aim*uMre*vMre*vMim+dt*aim*uMim*vMre*vMre-dt*aim*uMim*vMim*vMim;

uim[ind]=uoldim-dt*aim*uMre*vMre*vMre+dt*are*uMim*vMim*vMim+2.0*dt*aim*uMim*vMre*vMim-2.0*dt*are*uMre*vMre*vMim-dt*are*uMim*vMre*vMre+dt*aim*uMre*vMim*vMim;

vre[ind]=voldre+dt*are*uMre*vMre*vMre-dt*are*uMre*vMim*vMim-2.0*dt*are*uMim*vMre*vMim-2.0*dt*aim*uMre*vMre*vMim-dt*aim*uMim*vMre*vMre+dt*aim*uMim*vMim*vMim;

vim[ind]=voldim+dt*aim*uMre*vMre*vMre-dt*are*uMim*vMim*vMim-2.0*dt*aim*uMim*vMre*vMim+2.0*dt*are*uMre*vMre*vMim+dt*are*uMim*vMre*vMre-dt*aim*uMre*vMim*vMim;
							
chg=sqrt((ure[ind]-utempre)*(ure[ind]-utempre)+(uim[ind]-utempim)*(uim[ind]-utempim))+
    sqrt((vre[ind]-vtempre)*(vre[ind]-vtempre)+(vim[ind]-vtempim)*(vim[ind]-vtempim));
}   	
}
