#include "clFFT.h"
//#include <CL/cl.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#if defined(cl_amd_fp64) || defined(cl_khr_fp64)
    #if defined(cl_amd_fp64)
    #pragma OPENCL EXTENSION cl_amd_fp64 : enable
    #elif defined(cl_khr_fp64)
    #pragma OPENCL EXTENSION cl_khr_fp64 : enable
    #endif
// function declarations/definitions using double precision doubleing-point arithmetic
#endif
#define MAX_SOURCE_SIZE (0x100000)

int main(void) {
	//int	  reclen;
	int	  Nx;
	int 	  Ny;
	int 	  Nz;
	int	  N;
	int mm;
	double	  Tmax=0.0;
	double	  plotgap=0.0;
	double	  Lx,Ly,Lz;
	double	  dt=0.0;	
	double	  A=0.0;
	double	  B=0.0;
	double	  Du=0.0;
	double	  Dv=0.0;
	double	  errortol=0.0;
	double*	  uhigh[2],*vhigh[2];
	clock_t begin, end;
	double*		  time,*error;
	int*		  tries;// !#of times a timestep has been done
	double	  allerr[2];
	//! splitting coeffiecents
	double ahigh[4][2]={{0.0,0.0},
		{0.3243964040201712,0.1345862724908067},
		{0.3512071919596576,-0.2691725449816134},
		{0.3243964040201712,0.1345862724908067}};
	double  bhigh[4][2]={{0.1621982020100856,0.0672931362454034},
		{0.3378017979899144,-0.0672931362454034},
		{0.3378017979899144,-0.0672931362454034},
		{0.1621982020100856,0.0672931362454034}};

//openCL variables and initialize opencl
    cl_platform_id platform_id = NULL;
    cl_device_id device_id = NULL;
    cl_context context = NULL;
    cl_command_queue command_queue = NULL;
    cl_mem cl_uhigh[2] = {NULL,NULL};
    cl_mem cl_vhigh[2] = {NULL,NULL};
    cl_mem cl_uhat[2] = {NULL,NULL};
    cl_mem cl_vhat[2] = {NULL,NULL};
    cl_mem cl_uoldstep[2] = {NULL,NULL};
    cl_mem cl_voldstep[2] = {NULL,NULL};
    cl_mem cl_ulow[2] = {NULL,NULL};
    cl_mem cl_vlow[2] = {NULL,NULL};
    cl_mem cl_x = NULL;
    cl_mem cl_y = NULL;
    cl_mem cl_z = NULL;
    cl_mem cl_kx = NULL;
    cl_mem cl_ky = NULL;
    cl_mem cl_kz = NULL;
    cl_program p_grid = NULL,p_frequencies = NULL,p_initialdata = NULL, p_abserr=NULL,p_linearpart=NULL,p_nonlinearpart=NULL;//p_reduceerr=NULL,
    cl_kernel grid = NULL,frequencies = NULL,initialdata = NULL, abserr=NULL,linearpart=NULL,nonlinearpart=NULL;//reduceerr=NULL,
    cl_uint ret_num_devices;
    cl_uint ret_num_platforms;
    cl_int ret;
	ret = clGetPlatformIDs(1, &platform_id, &ret_num_platforms);
    	ret = clGetDeviceIDs(platform_id, CL_DEVICE_TYPE_DEFAULT, 1, &device_id, &ret_num_devices);
	context = clCreateContext( NULL, 1, &device_id, NULL, NULL, &ret);
	command_queue = clCreateCommandQueue(context, device_id, 0, &ret);
    	size_t source_size;
    	char *source_str;
//end opencl
	int  i,l,n,ierr;
	int status=0;	
//Read infutfile
	char	InputFileName[]="./INPUTFILEa";
	FILE*fp;
	fp=fopen(InputFileName,"r");
   	 if(!fp) {fprintf(stderr, "Failed to load IPUTFILEa.\n");exit(1);}	 
	ierr=fscanf(fp, "%d %d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", &Nx,&Ny,&Nz,&Tmax,&plotgap,&Lx,&Ly,&Lz,&dt,&Du,&Dv,&A,&B,&errortol);
	if(ierr!=14){fprintf(stderr, "INPUTFILE corrupted.\n");exit(1);}	
	fclose(fp);
	printf("NX %d\n",Nx); 
	printf("NY %d\n",Ny); 
	printf("NZ %d\n",Nz); 
	printf("Tmax %lf\n",Tmax);
	printf("plotgap %lf\n",plotgap);
	printf("Lx %lf\n",Lx);
	printf("Ly %lf\n",Ly);
	printf("Lz %lf\n",Lz);
	printf("dt %lf\n",dt);		
	printf("Du %lf\n",Du);
	printf("Dv %lf\n",Dv);
	printf("F %lf\n",A);
	printf("k %lf\n",B	);
	printf("errortol %lf\n",errortol);
	printf("Read inputfile\n");
	N=Nx*Ny*Nz;
	B=A+B;
//ALLocate the memory
	uhigh[0]=(double*) malloc(N*sizeof(double));
	vhigh[0]=(double*) malloc(N*sizeof(double));
	time=(double*) malloc(1000000*sizeof(double));
	error=(double*) malloc(1000000*sizeof(double));
	tries=(int*) malloc(1000000*sizeof(int));
//allocat gpu mem	
	cl_uhigh[0] = clCreateBuffer(context, CL_MEM_READ_WRITE, N * sizeof(double), NULL, &ret);
	cl_vhigh[0] = clCreateBuffer(context, CL_MEM_READ_WRITE, N * sizeof(double), NULL, &ret);
	cl_uhigh[1] = clCreateBuffer(context, CL_MEM_READ_WRITE, N * sizeof(double), NULL, &ret);
	cl_vhigh[1] = clCreateBuffer(context, CL_MEM_READ_WRITE, N * sizeof(double), NULL, &ret);
	cl_uoldstep[0] = clCreateBuffer(context, CL_MEM_READ_WRITE, N * sizeof(double), NULL, &ret);
	cl_voldstep[0] = clCreateBuffer(context, CL_MEM_READ_WRITE, N * sizeof(double), NULL, &ret);
	cl_uoldstep[1] = clCreateBuffer(context, CL_MEM_READ_WRITE, N * sizeof(double), NULL, &ret);
	cl_voldstep[1] = clCreateBuffer(context, CL_MEM_READ_WRITE, N * sizeof(double), NULL, &ret);
	cl_uhat[0] = clCreateBuffer(context, CL_MEM_READ_WRITE, N * sizeof(double), NULL, &ret);
	cl_vhat[0] = clCreateBuffer(context, CL_MEM_READ_WRITE, N * sizeof(double), NULL, &ret);
	cl_uhat[1] = clCreateBuffer(context, CL_MEM_READ_WRITE, N * sizeof(double), NULL, &ret);
	cl_vhat[1] = clCreateBuffer(context, CL_MEM_READ_WRITE, N * sizeof(double), NULL, &ret);
	cl_ulow[0] = clCreateBuffer(context, CL_MEM_READ_WRITE, N * sizeof(double), NULL, &ret);
	cl_vlow[0] = clCreateBuffer(context, CL_MEM_READ_WRITE, N * sizeof(double), NULL, &ret);
	cl_ulow[1] = clCreateBuffer(context, CL_MEM_READ_WRITE, N * sizeof(double), NULL, &ret);
	cl_vlow[1] = clCreateBuffer(context, CL_MEM_READ_WRITE, N * sizeof(double), NULL, &ret);
	printf("allocated space\n");

// FFT library realted declarations. 
	clfftPlanHandle planHandle;
	clfftDim dim = CLFFT_3D;
	size_t clLengths[3] = {Nx, Ny, Nz};
	// Setup clFFT. 
	clfftSetupData fftSetup;
	ret = clfftInitSetupData(&fftSetup);
	ret = clfftSetup(&fftSetup);
	// Create a default plan for a complex FFT. 
	ret = clfftCreateDefaultPlan(&planHandle, context, dim, clLengths);
	// Set plan parameters. 
	ret = clfftSetPlanPrecision(planHandle, CLFFT_DOUBLE);
	ret = clfftSetLayout(planHandle, CLFFT_COMPLEX_PLANAR, CLFFT_COMPLEX_PLANAR);
	ret = clfftSetResultLocation(planHandle, CLFFT_OUTOFPLACE);
	// Bake the plan. 
	ret = clfftBakePlan(planHandle, 1, &command_queue, NULL, NULL);
	// Create temporary buffer. 
	cl_mem tmpBufferu = 0;
	cl_mem tmpBufferv = 0;
	// Size of temp buffer. 
	size_t tmpBufferSize = 0;
	status = clfftGetTmpBufSize(planHandle, &tmpBufferSize);
	if ((status == 0) && (tmpBufferSize > 0)) {
		tmpBufferu = clCreateBuffer(context, CL_MEM_READ_WRITE, tmpBufferSize, NULL, &ret);
		tmpBufferv = clCreateBuffer(context, CL_MEM_READ_WRITE, tmpBufferSize, NULL, &ret);
		if (ret != CL_SUCCESS)
			printf("Error with tmpBuffer clCreateBuffer\n");
	}

//kernel grid
    	fp = fopen("./grid.c", "r");
    	if (!fp) {fprintf(stderr, "Failed to load grid.\n"); exit(1); }
    	source_str = (char *)malloc(MAX_SOURCE_SIZE);
   	source_size = fread( source_str, 1, MAX_SOURCE_SIZE, fp );
    	fclose( fp );
	
	p_grid = clCreateProgramWithSource(context, 1, (const char **)&source_str, (const size_t *)&source_size, &ret);
        ret = clBuildProgram(p_grid, 1, &device_id, NULL, NULL, NULL);
        grid = clCreateKernel(p_grid, "grid", &ret);
//get grid first x
	cl_x = clCreateBuffer(context, CL_MEM_READ_WRITE, Nx * sizeof(double), NULL, &ret);
        ret = clSetKernelArg(grid, 0, sizeof(cl_mem), (void *)&cl_x);
	ret = clSetKernelArg(grid, 1, sizeof(double),(void*)&Lx);
	ret = clSetKernelArg(grid, 2, sizeof(int),(void*)&Nx);
	size_t global_work_size_x[3] = {Nx, 0, 0};
        ret = clEnqueueNDRangeKernel(command_queue, grid, 1, NULL, global_work_size_x, 0, 0, NULL, NULL);
	ret = clFinish(command_queue);
//then y
	cl_y = clCreateBuffer(context, CL_MEM_READ_WRITE, Ny * sizeof(double), NULL, &ret);	
	ret = clSetKernelArg(grid, 0, sizeof(cl_mem), (void *)&cl_y);
	ret = clSetKernelArg(grid, 1, sizeof(double),(void*)&Ly);
	ret = clSetKernelArg(grid, 2, sizeof(int),(void*)&Ny);
	size_t global_work_size_y[3] = {Ny, 0, 0};
	ret = clEnqueueNDRangeKernel(command_queue, grid, 1, NULL, global_work_size_y, 0, 0, NULL, NULL);
	ret = clFinish(command_queue);


//last z
	cl_z = clCreateBuffer(context, CL_MEM_READ_WRITE, Nz * sizeof(double), NULL, &ret);
	ret = clSetKernelArg(grid, 0, sizeof(cl_mem), (void *)&cl_z);
	ret = clSetKernelArg(grid, 1, sizeof(double),(void*)&Lz);
	ret = clSetKernelArg(grid, 2, sizeof(int),(void*)&Nz);
	size_t global_work_size_z[3] = {Nz, 0, 0};
	ret = clEnqueueNDRangeKernel(command_queue, grid, 1, NULL, global_work_size_z, 0, 0, NULL, NULL);
	ret = clFinish(command_queue);
    	ret = clReleaseKernel(grid); ret = clReleaseProgram(p_grid);

//kernel initial data
    	fp = fopen("./initialdata.c", "r");
    	if (!fp) {fprintf(stderr, "Failed to load initialdata.\n"); exit(1); }
	free(source_str);    	
	source_str = (char *)malloc(MAX_SOURCE_SIZE);
   	source_size = fread( source_str, 1, MAX_SOURCE_SIZE, fp );
    	fclose( fp );

	p_initialdata = clCreateProgramWithSource(context, 1, (const char **)&source_str, (const size_t *)&source_size, &ret);
        ret = clBuildProgram(p_initialdata, 1, &device_id, NULL, NULL, NULL);
        initialdata = clCreateKernel(p_initialdata, "initialdata", &ret);
	size_t global_work_size[3] = {N, 0, 0};



//frequencies kernel

    	fp = fopen("./frequencies.c", "r");
    	if (!fp) {fprintf(stderr, "Failed to load frequencies.\n"); exit(1); }
	free(source_str);
    	source_str = (char *)malloc(MAX_SOURCE_SIZE);
   	source_size = fread( source_str, 1, MAX_SOURCE_SIZE, fp );
    	fclose( fp );
	
	p_frequencies = clCreateProgramWithSource(context, 1, (const char **)&source_str, (const size_t *)&source_size, &ret);
        ret = clBuildProgram(p_frequencies, 1, &device_id, NULL, NULL, NULL);
        frequencies = clCreateKernel(p_frequencies, "frequencies", &ret);
//get frequencies first x
	cl_kx = clCreateBuffer(context, CL_MEM_READ_WRITE, Nx * sizeof(double), NULL, &ret);
        ret = clSetKernelArg(frequencies, 0, sizeof(cl_mem), (void *)&cl_kx);
	ret = clSetKernelArg(frequencies, 1, sizeof(double),(void*)&Lx);
	ret = clSetKernelArg(frequencies, 2, sizeof(int),(void*)&Nx);
        ret = clEnqueueNDRangeKernel(command_queue, frequencies, 1, NULL, global_work_size_x, 0, 0, NULL, NULL);
	ret = clFinish(command_queue);
//then y
	cl_ky = clCreateBuffer(context, CL_MEM_READ_WRITE, Ny * sizeof(double), NULL, &ret);	
	ret = clSetKernelArg(frequencies, 0, sizeof(cl_mem), (void *)&cl_ky);
	ret = clSetKernelArg(frequencies, 1, sizeof(double),(void*)&Ly);
	ret = clSetKernelArg(frequencies, 2, sizeof(int),(void*)&Ny);
	ret = clEnqueueNDRangeKernel(command_queue, frequencies, 1, NULL, global_work_size_y, 0, 0, NULL, NULL);
	ret = clFinish(command_queue);
//last z
	cl_kz = clCreateBuffer(context, CL_MEM_READ_WRITE, Nz * sizeof(double), NULL, &ret);
	ret = clSetKernelArg(frequencies, 0, sizeof(cl_mem), (void *)&cl_kz);
	ret = clSetKernelArg(frequencies, 1, sizeof(double),(void*)&Lz);
	ret = clSetKernelArg(frequencies, 2, sizeof(int),(void*)&Nz);
	ret = clEnqueueNDRangeKernel(command_queue, frequencies, 1, NULL, global_work_size_z, 0, 0, NULL, NULL);
	ret = clFinish(command_queue);

	printf("Setup grid, fourier frequencies and initialcondition\n");
//load the rest of the kernels
//linearpart kernel
    	fp = fopen("./linearpart.c", "r");
    	if (!fp) {fprintf(stderr, "Failed to load linearpart.\n"); exit(1); }
	free(source_str);    	
	source_str = (char *)malloc(MAX_SOURCE_SIZE);
   	source_size = fread( source_str, 1, MAX_SOURCE_SIZE, fp );
    	fclose( fp );

	p_linearpart = clCreateProgramWithSource(context, 1, (const char **)&source_str, (const size_t *)&source_size, &ret);
        ret = clBuildProgram(p_linearpart, 1, &device_id, NULL, NULL, NULL);
        linearpart = clCreateKernel(p_linearpart, "linearpart", &ret);

//kernel nonlinear
    	fp = fopen("./nonlinearpart.c", "r");
    	if (!fp) {fprintf(stderr, "Failed to load nonlinearpart.\n"); exit(1); }
	free(source_str);    	
	source_str = (char *)malloc(MAX_SOURCE_SIZE);
   	source_size = fread( source_str, 1, MAX_SOURCE_SIZE, fp );
    	fclose( fp );

	p_nonlinearpart = clCreateProgramWithSource(context, 1, (const char **)&source_str, (const size_t *)&source_size, &ret);
        ret = clBuildProgram(p_nonlinearpart, 1, &device_id, NULL, NULL, NULL);
        nonlinearpart = clCreateKernel(p_nonlinearpart, "nonlinearpart", &ret);

//kernel abserr
	fp = fopen("./abserr.c", "r");
    	if (!fp) {fprintf(stderr, "Failed to load abserr.\n"); exit(1); }
	free(source_str);    	
	source_str = (char *)malloc(MAX_SOURCE_SIZE);
   	source_size = fread( source_str, 1, MAX_SOURCE_SIZE, fp );
    	fclose( fp );

	p_abserr = clCreateProgramWithSource(context, 1, (const char **)&source_str, (const size_t *)&source_size, &ret);
        ret = clBuildProgram(p_abserr, 1, &device_id, NULL, NULL, NULL);
        abserr = clCreateKernel(p_abserr, "abserr", &ret);

//initialise oldstep
for(mm=0;mm<=12;mm++){
//get initial data
        ret = clSetKernelArg(initialdata, 0, sizeof(cl_mem),(void *)&cl_uhigh[0]);
	ret = clSetKernelArg(initialdata, 1, sizeof(cl_mem),(void* )&cl_vhigh[0]);
        ret = clSetKernelArg(initialdata, 2, sizeof(cl_mem),(void *)&cl_uhigh[1]);
	ret = clSetKernelArg(initialdata, 3, sizeof(cl_mem),(void* )&cl_vhigh[1]);
	ret = clSetKernelArg(initialdata, 4, sizeof(cl_mem),(void* )&cl_x);
	ret = clSetKernelArg(initialdata, 5, sizeof(cl_mem),(void* )&cl_y);
	ret = clSetKernelArg(initialdata, 6, sizeof(cl_mem),(void* )&cl_z);
	ret = clSetKernelArg(initialdata, 7, sizeof(int),(void* )&Nx);
	ret = clSetKernelArg(initialdata, 8, sizeof(int),(void* )&Ny);
	ret = clSetKernelArg(initialdata, 9, sizeof(int),(void* )&Nz);

        ret = clEnqueueNDRangeKernel(command_queue, initialdata, 1, NULL, global_work_size, 0, 0, NULL, NULL);
	ret = clFinish(command_queue);

	n=0;
	time[n]=0.0;
	printf("Got initial data, starting timestepping\n");
  	begin = clock();
	while(time[n]<Tmax){
	n=n+1;
	time[n]=time[n-1]+dt;
//linear
	ret = clfftEnqueueTransform(planHandle, CLFFT_FORWARD, 1, &command_queue, 0, NULL, NULL,cl_uhigh, cl_uhat, tmpBufferu);
	ret = clfftEnqueueTransform(planHandle, CLFFT_FORWARD, 1, &command_queue, 0, NULL, NULL,cl_vhigh, cl_vhat, tmpBufferv);
	ret = clFinish(command_queue);

        ret = clSetKernelArg(linearpart, 0, sizeof(cl_mem),(void *)&cl_uhat[0]);
        ret = clSetKernelArg(linearpart, 1, sizeof(cl_mem),(void *)&cl_uhat[1]);
        ret = clSetKernelArg(linearpart, 2, sizeof(cl_mem),(void *)&cl_vhat[0]);
        ret = clSetKernelArg(linearpart, 3, sizeof(cl_mem),(void *)&cl_vhat[1]);
	ret = clSetKernelArg(linearpart, 4, sizeof(cl_mem),(void* )&cl_kx);
	ret = clSetKernelArg(linearpart, 5, sizeof(cl_mem),(void* )&cl_ky);
	ret = clSetKernelArg(linearpart, 6, sizeof(cl_mem),(void* )&cl_kz);
	ret = clSetKernelArg(linearpart, 7, sizeof(double),(void* )&dt);
	ret = clSetKernelArg(linearpart, 8, sizeof(double),(void* )&Du);
	ret = clSetKernelArg(linearpart, 9, sizeof(double),(void* )&Dv);
	ret = clSetKernelArg(linearpart, 10, sizeof(double),(void* )&A);
	ret = clSetKernelArg(linearpart, 11, sizeof(double),(void* )&B);
	ret = clSetKernelArg(linearpart, 12, sizeof(double),(void* )&bhigh[0][0]);
	ret = clSetKernelArg(linearpart, 13, sizeof(double),(void* )&bhigh[0][1]);
	ret = clSetKernelArg(linearpart, 14, sizeof(int),(void* )&Nx);
	ret = clSetKernelArg(linearpart, 15, sizeof(int),(void* )&Ny);
	ret = clSetKernelArg(linearpart, 16, sizeof(int),(void* )&Nz);
        ret = clEnqueueNDRangeKernel(command_queue, linearpart, 1, NULL, global_work_size, 0, 0, NULL, NULL);
	ret = clFinish(command_queue);

	ret = clfftEnqueueTransform(planHandle, CLFFT_BACKWARD, 1, &command_queue, 0, NULL, NULL,cl_uhat, cl_uhigh, tmpBufferu);
	ret = clfftEnqueueTransform(planHandle, CLFFT_BACKWARD, 1, &command_queue, 0, NULL, NULL,cl_vhat, cl_vhigh, tmpBufferv);
	ret = clFinish(command_queue);
//!a(1),b(1) done

for(l=1;l<4;l++){
//nonlinearpart
        ret = clSetKernelArg(nonlinearpart, 0, sizeof(cl_mem),(void *)&cl_uhigh[0]);
        ret = clSetKernelArg(nonlinearpart, 1, sizeof(cl_mem),(void *)&cl_uhigh[1]);
	ret = clSetKernelArg(nonlinearpart, 2, sizeof(cl_mem),(void* )&cl_vhigh[0]);
	ret = clSetKernelArg(nonlinearpart, 3, sizeof(cl_mem),(void* )&cl_vhigh[1]);
	ret = clSetKernelArg(nonlinearpart, 4, sizeof(double),(void* )&dt);
	ret = clSetKernelArg(nonlinearpart, 5, sizeof(double),(void* )&ahigh[l][0]);
	ret = clSetKernelArg(nonlinearpart, 6, sizeof(double),(void* )&ahigh[l][1]);
        ret = clEnqueueNDRangeKernel(command_queue, nonlinearpart, 1, NULL, global_work_size, 0, 0, NULL, NULL);
	ret = clFinish(command_queue);
//linear
	ret = clfftEnqueueTransform(planHandle, CLFFT_FORWARD, 1, &command_queue, 0, NULL, NULL,cl_uhigh, cl_uhat, tmpBufferu);
	ret = clfftEnqueueTransform(planHandle, CLFFT_FORWARD, 1, &command_queue, 0, NULL, NULL,cl_vhigh, cl_vhat, tmpBufferv);
	ret = clFinish(command_queue);

        ret = clSetKernelArg(linearpart, 0, sizeof(cl_mem),(void *)&cl_uhat[0]);
        ret = clSetKernelArg(linearpart, 1, sizeof(cl_mem),(void *)&cl_uhat[1]);
        ret = clSetKernelArg(linearpart, 2, sizeof(cl_mem),(void *)&cl_vhat[0]);
        ret = clSetKernelArg(linearpart, 3, sizeof(cl_mem),(void *)&cl_vhat[1]);
	ret = clSetKernelArg(linearpart, 4, sizeof(cl_mem),(void* )&cl_kx);
	ret = clSetKernelArg(linearpart, 5, sizeof(cl_mem),(void* )&cl_ky);
	ret = clSetKernelArg(linearpart, 6, sizeof(cl_mem),(void* )&cl_kz);
	ret = clSetKernelArg(linearpart, 7, sizeof(double),(void* )&dt);
	ret = clSetKernelArg(linearpart, 8, sizeof(double),(void* )&Du);
	ret = clSetKernelArg(linearpart, 9, sizeof(double),(void* )&Dv);
	ret = clSetKernelArg(linearpart, 10, sizeof(double),(void* )&A);
	ret = clSetKernelArg(linearpart, 11, sizeof(double),(void* )&B);
	ret = clSetKernelArg(linearpart, 12, sizeof(double),(void* )&bhigh[l][0]);
	ret = clSetKernelArg(linearpart, 13, sizeof(double),(void* )&bhigh[l][1]);
	ret = clSetKernelArg(linearpart, 14, sizeof(int),(void* )&Nx);
	ret = clSetKernelArg(linearpart, 15, sizeof(int),(void* )&Ny);
	ret = clSetKernelArg(linearpart, 16, sizeof(int),(void* )&Nz);
        ret = clEnqueueNDRangeKernel(command_queue, linearpart, 1, NULL, global_work_size, 0, 0, NULL, NULL);
	ret = clFinish(command_queue);

	ret = clfftEnqueueTransform(planHandle, CLFFT_BACKWARD, 1, &command_queue, 0, NULL, NULL,cl_uhat, cl_uhigh, tmpBufferu);
	ret = clfftEnqueueTransform(planHandle, CLFFT_BACKWARD, 1, &command_queue, 0, NULL, NULL,cl_vhat, cl_vhigh, tmpBufferv);
	ret = clFinish(command_queue);	}
//high done
	}
		ret=clEnqueueCopyBuffer(command_queue,cl_uoldstep[0],cl_ulow[0],0,0,N*sizeof(double),0,NULL,NULL);
		ret=clEnqueueCopyBuffer(command_queue,cl_voldstep[0],cl_vlow[0],0,0,N*sizeof(double),0,NULL,NULL);
		ret=clEnqueueCopyBuffer(command_queue,cl_uoldstep[1],cl_ulow[1],0,0,N*sizeof(double),0,NULL,NULL);
		ret=clEnqueueCopyBuffer(command_queue,cl_voldstep[1],cl_vlow[1],0,0,N*sizeof(double),0,NULL,NULL);
		ret = clFinish(command_queue);
	//uerror
        ret = clSetKernelArg(abserr, 0, sizeof(cl_mem),(void *)&cl_uhigh[0]);
        ret = clSetKernelArg(abserr, 1, sizeof(cl_mem),(void *)&cl_uhigh[1]);
	ret = clSetKernelArg(abserr, 2, sizeof(cl_mem),(void* )&cl_ulow[0]);
	ret = clSetKernelArg(abserr, 3, sizeof(cl_mem),(void* )&cl_ulow[1]);
	ret = clEnqueueNDRangeKernel(command_queue, abserr, 1, NULL, global_work_size, 0, 0, NULL, NULL);
	ret = clFinish(command_queue);
        ret = clEnqueueReadBuffer(command_queue, cl_ulow[0], CL_TRUE, 0, N * sizeof(double), uhigh[0], 0, NULL, NULL);
	ret = clFinish(command_queue);
 	allerr[0]=0.0;
	int zahler=0;
	for(i=0;i<N;i++){
		if(uhigh[0][i]>allerr[0]){allerr[0]=uhigh[0][i];zahler++;}
	}
//verror
	ret = clSetKernelArg(abserr, 0, sizeof(cl_mem),(void *)&cl_vhigh[0]);
        ret = clSetKernelArg(abserr, 1, sizeof(cl_mem),(void *)&cl_vhigh[1]);
	ret = clSetKernelArg(abserr, 2, sizeof(cl_mem),(void* )&cl_vlow[0]);
	ret = clSetKernelArg(abserr, 3, sizeof(cl_mem),(void* )&cl_vlow[1]);
	ret = clEnqueueNDRangeKernel(command_queue, abserr, 1, NULL, global_work_size, 0, 0, NULL, NULL);
	ret = clFinish(command_queue);
        ret = clEnqueueReadBuffer(command_queue, cl_vlow[0], CL_TRUE, 0, N * sizeof(double), uhigh[0], 0, NULL, NULL);
	ret = clFinish(command_queue);
 	allerr[1]=0.0;
 	zahler=0;
	for(i=0;i<N;i++){
//printf("%lf ",uhigh[0][i]);//-uhigh[0][i]);
	if(uhigh[0][i]>allerr[1]){allerr[1]=uhigh[0][i];zahler++;}
	}

	error[mm]=allerr[0]+allerr[1];
//prepare for next time step
		ret=clEnqueueCopyBuffer(command_queue,cl_uhigh[0],cl_uoldstep[0],0,0,N*sizeof(double),0,NULL,NULL);
		ret=clEnqueueCopyBuffer(command_queue,cl_vhigh[0],cl_voldstep[0],0,0,N*sizeof(double),0,NULL,NULL);
		ret=clEnqueueCopyBuffer(command_queue,cl_uhigh[1],cl_uoldstep[1],0,0,N*sizeof(double),0,NULL,NULL);
		ret=clEnqueueCopyBuffer(command_queue,cl_vhigh[1],cl_voldstep[1],0,0,N*sizeof(double),0,NULL,NULL);
		ret = clFinish(command_queue);
	end = clock();
	printf("Finished time stepping\n");
	printf("%d,Programm took %lf for execution, error%e\n",mm,(double)(end-begin)/CLOCKS_PER_SEC,error[mm]);
// Save error
	fp=fopen("./data/edata43.dat","w");
    	if (!fp) {fprintf(stderr, "Failed to write edata.dat.\n"); exit(1); }
	for(i=0;i<=mm;i++){fprintf(fp,"%e\n",error[i]);}
    	fclose( fp );
	dt=dt*0.5;
}
//release memory
	clReleaseMemObject(cl_uhigh[0]);
	clReleaseMemObject(cl_uhigh[1]);
	clReleaseMemObject(cl_vhigh[0]);
	clReleaseMemObject(cl_vhigh[1]);
	clReleaseMemObject(cl_uoldstep[0]);
	clReleaseMemObject(cl_uoldstep[1]);
	clReleaseMemObject(cl_voldstep[0]);
	clReleaseMemObject(cl_voldstep[1]);
	clReleaseMemObject(cl_uhat[0]);
	clReleaseMemObject(cl_uhat[1]);
	clReleaseMemObject(cl_vhat[0]);
	clReleaseMemObject(cl_vhat[1]);
	clReleaseMemObject(cl_kx);
	clReleaseMemObject(cl_ky);
	clReleaseMemObject(cl_kz);
	ret = clReleaseMemObject(cl_x);
	ret = clReleaseMemObject(cl_y);
	ret = clReleaseMemObject(cl_z);
    	ret = clReleaseKernel(frequencies); ret = clReleaseProgram(p_frequencies);
    	ret = clReleaseKernel(linearpart); ret = clReleaseProgram(p_linearpart);
    	ret = clReleaseKernel(nonlinearpart); ret = clReleaseProgram(p_nonlinearpart);
    	ret = clReleaseKernel(abserr); ret = clReleaseProgram(p_abserr);
	ret = clReleaseKernel(initialdata); ret = clReleaseProgram(p_initialdata);
    	//ret = clReleaseKernel(reduceerr); ret = clReleaseProgram(p_reduceerr);
	free(uhigh[0]);
	free(vhigh[0]);
	free(time);
	free(error);
	free(tries);
//release clFFT
	clReleaseMemObject(tmpBufferu);
	clReleaseMemObject(tmpBufferv);
	/* Release the plan. */
	ret = clfftDestroyPlan(&planHandle);
	/* Release clFFT library. */
	clfftTeardown();
//release OpenCL
	ret = clReleaseCommandQueue(command_queue);
     	ret = clReleaseContext(context);	
	printf("Program execution complete\n");
	return 0;
}
