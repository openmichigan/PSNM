#include "clFFT.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include "gs_functions.c"

#if defined(cl_amd_fp64) || defined(cl_khr_fp64)
    #if defined(cl_amd_fp64)
    #pragma OPENCL EXTENSION cl_amd_fp64 : enable
    #elif defined(cl_khr_fp64)
    #pragma OPENCL EXTENSION cl_khr_fp64 : enable
    #endif
    // function declarations/definitions using double precision floating-point arithmetic
#endif

int main(void) {
//time meassuring
  	struct timeval tvs;

//variables
	int 	Nx=1024;
	int		Ny=1024;
	int 	plotnum=0;
	int	  	Tmax=2;
	int 	plottime=0;
	int	  	plotgap=1;
	double	Lx=1.0;
	double 	Ly=1.0;
	double	dt=0.0;	
	double	A=0.0;
	double	B=0.0;
	double	Du=0.0;
	double	Dv=0.0;
//splitting coefficients
	double	a=0.5;	
	double 	b=0.5;
	double 	c=1.0;
//loop counters	
	int i=0;
	int j=0;
	int n=0;

	double*umax=NULL;
	double*vmax=NULL;
	parainit(&Nx,&Ny,&Tmax,&plotgap,&Lx,&Ly,&dt,&Du,&Dv,&A,&B);
	plottime=plotgap;
	vmax=(double*)malloc((Tmax/plotgap+1)*sizeof(double));
	umax=(double*)malloc((Tmax/plotgap+1)*sizeof(double));
//openCL variables
    cl_platform_id *platform_id = NULL;
    cl_kernel frequencies = NULL, initialdata = NULL, linearpart=NULL;
	cl_kernel nonlinearpart_a=NULL, nonlinearpart_b=NULL;
    cl_int ret;
    cl_uint num_platforms;
// Detect how many platforms there are.
	ret = clGetPlatformIDs(0, NULL, &num_platforms);
// Allocate enough space for the number of platforms.
	platform_id = (cl_platform_id*) malloc(num_platforms*sizeof(cl_platform_id));
// Store the platforms
	ret = clGetPlatformIDs(num_platforms, platform_id, NULL);
	printf("Found %d platform(s)!\n",num_platforms);
    cl_uint *num_devices;
	num_devices=(cl_uint*) malloc(num_platforms*sizeof(cl_uint));
    cl_device_id **device_id = NULL;
	device_id =(cl_device_id**) malloc(num_platforms*sizeof(cl_device_id*));
// Detect number of devices in the platforms
	for(i=0;i<num_platforms;i++){
		char buf[65536];
		size_t size;
		ret = clGetPlatformInfo(platform_id[i],CL_PLATFORM_VERSION,sizeof(buf),buf,&size);
		printf("%s\n",buf);
		ret = clGetDeviceIDs(platform_id[i],CL_DEVICE_TYPE_ALL,0,NULL,num_devices);
		printf("Found %d device(s) on platform %d!\n", num_devices[i],i);
		ret = clGetPlatformInfo(platform_id[i],CL_PLATFORM_NAME,sizeof(buf),buf,&size);
		printf("%s ",buf);
// Store numDevices from platform
		device_id[i]=(cl_device_id*) malloc(num_devices[i]*sizeof(device_id));
		ret = clGetDeviceIDs(platform_id[i],CL_DEVICE_TYPE_ALL,num_devices[i],device_id[i],NULL);
		for(j=0;j<num_devices[i];j++){
			ret = clGetDeviceInfo(device_id[i][j],CL_DEVICE_NAME,sizeof(buf),buf,&size);
			printf("%s (%d,%d)\n",buf,i,j);
		}
	}
//create context and command_queue
    cl_context context = NULL;
   	cl_command_queue command_queue = NULL;
//Which platform and device do i choose?
	int	chooseplatform=0;
	int	choosedevice=0;	  
	printf("Choose platform %d and device %d!\n",chooseplatform,choosedevice);
	context = clCreateContext( NULL, num_devices[chooseplatform], device_id[chooseplatform], NULL, NULL, &ret);
	if(ret!=CL_SUCCESS){printf("createContext ret:%d\n",ret); exit(1); }
	command_queue = clCreateCommandQueue(context, device_id[chooseplatform][choosedevice], 0, &ret);
	if(ret!=CL_SUCCESS){printf("createCommandQueue ret:%d\n",ret); exit(1); }

//OpenCL arrays
    cl_mem cl_u = NULL,cl_v = NULL;
   	cl_mem cl_uhat = NULL, cl_vhat = NULL;
    cl_mem cl_kx = NULL, cl_ky = NULL;

//FFT
	clfftPlanHandle planHandle;
    cl_mem tmpBuffer = NULL;
	fftinit(&planHandle,&context, &command_queue, &tmpBuffer, Nx, Ny);

//allocate gpu memory/
	cl_u=clCreateBuffer(context, CL_MEM_READ_WRITE, 2*Nx* Ny* sizeof(double), NULL, &ret);
	cl_v=clCreateBuffer(context, CL_MEM_READ_WRITE, 2*Nx* Ny* sizeof(double), NULL, &ret);
	cl_uhat=clCreateBuffer(context, CL_MEM_READ_WRITE, 2*Nx * Ny* sizeof(double), NULL, &ret);
	cl_vhat=clCreateBuffer(context, CL_MEM_READ_WRITE, 2*Nx * Ny* sizeof(double), NULL, &ret);
	cl_kx = clCreateBuffer(context, CL_MEM_READ_WRITE, Nx * sizeof(double), NULL, &ret);
	cl_ky = clCreateBuffer(context, CL_MEM_READ_WRITE, Ny * sizeof(double), NULL, &ret);

	printf("allocated space\n");
//load the kernels
	loadKernel(&frequencies,&context,&device_id[chooseplatform][choosedevice],"frequencies");
	loadKernel(&initialdata,&context,&device_id[chooseplatform][choosedevice],"initialdata"); 
	loadKernel(&linearpart,&context,&device_id[chooseplatform][choosedevice],"linearpart"); 
	loadKernel(&nonlinearpart_a,&context,&device_id[chooseplatform][choosedevice],"nonlinearpart_a"); 
	loadKernel(&nonlinearpart_b,&context,&device_id[chooseplatform][choosedevice],"nonlinearpart_b"); 

	size_t global_work_size[1] = {Nx*Ny};
	size_t global_work_size_X[1] = {Nx};
	size_t global_work_size_Y[1] = {Ny};
//frequencies
    ret = clSetKernelArg(frequencies, 0, sizeof(cl_mem),(void *)&cl_kx);
	ret = clSetKernelArg(frequencies, 1, sizeof(double),(void* )&Lx);
	ret = clSetKernelArg(frequencies, 2, sizeof(int),(void* )&Nx);
    ret = clEnqueueNDRangeKernel(command_queue, frequencies, 1, NULL, global_work_size_X, NULL, 0, NULL, NULL);
	ret = clFinish(command_queue);
    ret = clSetKernelArg(frequencies, 0, sizeof(cl_mem),(void *)&cl_ky);
	ret = clSetKernelArg(frequencies, 1, sizeof(double),(void* )&Ly);
	ret = clSetKernelArg(frequencies, 2, sizeof(int),(void* )&Ny);
    ret = clEnqueueNDRangeKernel(command_queue, frequencies, 1, NULL, global_work_size_Y, NULL, 0, NULL, NULL);
	ret = clFinish(command_queue);
//printCL(&cl_kx,&command_queue,Nx,1);
//printCL(&cl_ky,&command_queue,1,Ny);
//inintial data
    ret = clSetKernelArg(initialdata, 0, sizeof(cl_mem),(void *)&cl_u);
	ret = clSetKernelArg(initialdata, 1, sizeof(cl_mem),(void* )&cl_v);
	ret = clSetKernelArg(initialdata, 2, sizeof(int),(void* )&Nx);
	ret = clSetKernelArg(initialdata, 3, sizeof(int),(void* )&Ny);
	ret = clSetKernelArg(initialdata, 4, sizeof(double),(void* )&Lx);
	ret = clSetKernelArg(initialdata, 5, sizeof(double),(void* )&Ly);
    ret = clEnqueueNDRangeKernel(command_queue, initialdata, 1, NULL, global_work_size, NULL, 0, NULL, NULL);
	ret = clFinish(command_queue);
//make output
    writedata_C(&cl_u, &command_queue,Nx,Ny,plotnum,"u");
    writedata_C(&cl_v, &command_queue,Nx,Ny,plotnum,"v");
    umax[plotnum]=writeimage(&cl_u, &command_queue,Nx,Ny,plotnum,"u");
    vmax[plotnum]=writeimage(&cl_v, &command_queue,Nx,Ny,plotnum,"v");
	printf("Got initial data, starting timestepping\n");
	mtime_s(&tvs);

	for(n=0;n<=Tmax;n++){
//nonlinearpart_a
    ret = clSetKernelArg(nonlinearpart_a, 0, sizeof(cl_mem),(void *)&cl_u);
	ret = clSetKernelArg(nonlinearpart_a, 1, sizeof(cl_mem),(void* )&cl_v);
	ret = clSetKernelArg(nonlinearpart_a, 2, sizeof(double),(void* )&A);
	ret = clSetKernelArg(nonlinearpart_a, 3, sizeof(double),(void* )&dt);
	ret = clSetKernelArg(nonlinearpart_a, 4, sizeof(double),(void* )&a);
    ret = clEnqueueNDRangeKernel(command_queue, nonlinearpart_a, 1, NULL, global_work_size, NULL, 0, NULL, NULL);
	ret = clFinish(command_queue);	

//nonlinearpart_b
    ret = clSetKernelArg(nonlinearpart_b, 0, sizeof(cl_mem),(void *)&cl_u);
	ret = clSetKernelArg(nonlinearpart_b, 1, sizeof(cl_mem),(void* )&cl_v);
	ret = clSetKernelArg(nonlinearpart_b, 2, sizeof(double),(void* )&A);
	ret = clSetKernelArg(nonlinearpart_b, 3, sizeof(double),(void* )&dt);
	ret = clSetKernelArg(nonlinearpart_b, 4, sizeof(double),(void* )&b);
    ret = clEnqueueNDRangeKernel(command_queue, nonlinearpart_b, 1, NULL, global_work_size, NULL, 0, NULL, NULL);
	ret = clFinish(command_queue);

//linear
	fft2dfor(&cl_u, &cl_uhat,&planHandle,&command_queue,&tmpBuffer);
	fft2dfor(&cl_v, &cl_vhat,&planHandle,&command_queue,&tmpBuffer);
//printf("A%f,B%f\n",A,B);
    ret = clSetKernelArg(linearpart, 0, sizeof(cl_mem),(void *)&cl_uhat);
    ret = clSetKernelArg(linearpart, 1, sizeof(cl_mem),(void *)&cl_vhat);
	ret = clSetKernelArg(linearpart, 2, sizeof(cl_mem),(void* )&cl_kx);
	ret = clSetKernelArg(linearpart, 3, sizeof(cl_mem),(void* )&cl_ky);
	ret = clSetKernelArg(linearpart, 4, sizeof(double),(void* )&Du);
	ret = clSetKernelArg(linearpart, 5, sizeof(double),(void* )&Dv);
	ret = clSetKernelArg(linearpart, 6, sizeof(double),(void* )&A);
	ret = clSetKernelArg(linearpart, 7, sizeof(double),(void* )&B);
	ret = clSetKernelArg(linearpart, 8, sizeof(double),(void* )&dt);
	ret = clSetKernelArg(linearpart, 9, sizeof(double),(void* )&c);
	ret = clSetKernelArg(linearpart, 10, sizeof(int),(void* )&Nx);
	ret = clSetKernelArg(linearpart, 11, sizeof(int),(void* )&Ny);
    ret = clEnqueueNDRangeKernel(command_queue, linearpart, 1, NULL, global_work_size, NULL, 0, NULL, NULL);
	ret = clFinish(command_queue);

	fft2dback(&cl_u, &cl_uhat,&planHandle,&command_queue,&tmpBuffer);
  	fft2dback(&cl_v, &cl_vhat,&planHandle,&command_queue,&tmpBuffer);

//nonlinearpart_b
    ret = clSetKernelArg(nonlinearpart_b, 0, sizeof(cl_mem),(void *)&cl_u);
	ret = clSetKernelArg(nonlinearpart_b, 1, sizeof(cl_mem),(void* )&cl_v);
	ret = clSetKernelArg(nonlinearpart_b, 2, sizeof(double),(void* )&A);
	ret = clSetKernelArg(nonlinearpart_b, 3, sizeof(double),(void* )&dt);
	ret = clSetKernelArg(nonlinearpart_b, 4, sizeof(double),(void* )&b);
    ret = clEnqueueNDRangeKernel(command_queue, nonlinearpart_b, 1, NULL, global_work_size, NULL, 0, NULL, NULL);
	ret = clFinish(command_queue);		
//nonlinearpart_a
    ret = clSetKernelArg(nonlinearpart_a, 0, sizeof(cl_mem),(void *)&cl_u);
	ret = clSetKernelArg(nonlinearpart_a, 1, sizeof(cl_mem),(void* )&cl_v);
	ret = clSetKernelArg(nonlinearpart_a, 2, sizeof(double),(void* )&A);
	ret = clSetKernelArg(nonlinearpart_a, 3, sizeof(double),(void* )&dt);
	ret = clSetKernelArg(nonlinearpart_a, 4, sizeof(double),(void* )&a);
    ret = clEnqueueNDRangeKernel(command_queue, nonlinearpart_a, 1, NULL, global_work_size, NULL, 0, NULL, NULL);
	ret = clFinish(command_queue);	
// done
	if(n==plottime){
		printf("time:%f, step:%d,%d,umax:%f,vmax:%f\n",n*dt,n,plotnum,umax[plotnum],vmax[plotnum]);
		plottime=plottime+plotgap;
		plotnum=plotnum+1;
   	 	writedata_C(&cl_u, &command_queue,Nx,Ny,plotnum,"u");
    	writedata_C(&cl_v, &command_queue,Nx,Ny,plotnum,"v");
        umax[plotnum]=writeimage(&cl_u, &command_queue,Nx,Ny,plotnum,"u");
        vmax[plotnum]=writeimage(&cl_v, &command_queue,Nx,Ny,plotnum,"v");
	}
}//end timestepping

	printf("Finished time stepping\n");
	mtime_e(&tvs,"Programm took:");
	writearray(umax,(Tmax/plotgap)+1,"u");
	writearray(vmax,(Tmax/plotgap)+1,"v");
	free(umax);
	free(vmax);	

	clReleaseMemObject(cl_u);
	clReleaseMemObject(cl_v);
	clReleaseMemObject(cl_uhat);
	clReleaseMemObject(cl_vhat);
	clReleaseMemObject(cl_kx);
	clReleaseMemObject(cl_ky);

    ret = clReleaseKernel(initialdata); 
    ret = clReleaseKernel(frequencies); 
    ret = clReleaseKernel(linearpart); 
    ret = clReleaseKernel(nonlinearpart_a);
    ret = clReleaseKernel(nonlinearpart_b);

	fftdestroy(&planHandle, &tmpBuffer);

	ret = clReleaseCommandQueue(command_queue);
    ret = clReleaseContext(context);

	for(i=0;i<num_platforms;i++){free(device_id[i]);}
	free(device_id);
	free(platform_id);
	free(num_devices);
	printf("Program execution complete\n");

	return 0;
}
