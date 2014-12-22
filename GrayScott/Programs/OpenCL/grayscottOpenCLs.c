#include "clFFT.h"
//#include <CL/cl.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/time.h>

#if defined(cl_amd_fp64) || defined(cl_khr_fp64)
    #if defined(cl_amd_fp64)
    #pragma OPENCL EXTENSION cl_amd_fp64 : enable
    #elif defined(cl_khr_fp64)
    #pragma OPENCL EXTENSION cl_khr_fp64 : enable
    #endif
    // function declarations/definitions using double precision floating-point arithmetic
#endif

#define MAX_SOURCE_SIZE (0x100000)

int main(void) {
//time meassuring
  	struct timeval tvs;
  	struct timeval tve;
    	double elapsedTime;
//variables
	int	  Nx;
	int 	  Ny;
	int 	  Nz;
	int	  N;
	int 	  plotnum=0;
	int	  Tmax=0;
	int 	  plottime=0;
	int	  plotgap=0;
	double	  Lx,Ly,Lz;
	double	  dt=0.0;	
	double	  A=0.0;
	double	  B=0.0;
	double	  Du=0.0;
	double	  Dv=0.0;
	double	  a[2]={1.0,0.0};	
	double 	  b[2]={0.5,0.0};
	double*	  x,*y,*z ;
	double*	  u[2],*v[2];
//openCL variables
    cl_platform_id platform_id = NULL;
    cl_device_id device_id = NULL;
    cl_context context = NULL;
    cl_command_queue command_queue = NULL;
    cl_mem cl_u[2] = {NULL,NULL};
    cl_mem cl_v[2] = {NULL,NULL};
    cl_mem cl_uhat[2] = {NULL,NULL};
    cl_mem cl_vhat[2] = {NULL,NULL};
    cl_mem cl_x = NULL;
    cl_mem cl_y = NULL;
    cl_mem cl_z = NULL;
    cl_mem cl_kx = NULL;
    cl_mem cl_ky = NULL;
    cl_mem cl_kz = NULL;
    cl_program p_grid = NULL,p_frequencies = NULL,p_initialdata = NULL,p_linearpart=NULL,p_nonlinearpart=NULL;
    cl_kernel grid = NULL,frequencies = NULL,initialdata = NULL,linearpart=NULL,nonlinearpart=NULL;
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
	int  i,n;
int status=0;	
 	//int  start, finish, count_rate, ind, numthreads
	char	nameconfig[100]="";
//Read infutfile
	char	InputFileName[]="./INPUTFILE";
	FILE*fp;
	fp=fopen(InputFileName,"r");
   	 if(!fp) {fprintf(stderr, "Failed to load IPUTFILE.\n");exit(1);}	 
	int ierr=fscanf(fp, "%d %d %d %d %d %lf %lf %lf %lf %lf %lf %lf %lf", &Nx,&Ny,&Nz,&Tmax,&plotgap,&Lx,&Ly,&Lz,&dt,&Du,&Dv,&A,&B);
	if(ierr!=13){fprintf(stderr, "INPUTFILE corrupted.\n");exit(1);}	
	fclose(fp);
	printf("NX %d\n",Nx); 
	printf("NY %d\n",Ny); 
	printf("NZ %d\n",Nz); 
	printf("Tmax %d\n",Tmax);
	printf("plotgap %d\n",plotgap);
	printf("Lx %lf\n",Lx);
	printf("Ly %lf\n",Ly);
	printf("Lz %lf\n",Lz);
	printf("dt %lf\n",dt);		
	printf("Du %lf\n",Du);
	printf("Dv %lf\n",Dv);
	printf("F %lf\n",A);
	printf("k %lf\n",B);
	printf("Read inputfile\n");
	N=Nx*Ny*Nz;
	plottime=plotgap;
	B=A+B;
//ALLocate the memory
	u[0]=(double*) malloc(N*sizeof(double));
	v[0]=(double*) malloc(N*sizeof(double));
	x=(double*) malloc(Nx*sizeof(double));
	y=(double*) malloc(Ny*sizeof(double));
	z=(double*) malloc(Nz*sizeof(double));

//allocate gpu mem
	cl_u[0] = clCreateBuffer(context, CL_MEM_READ_WRITE, N * sizeof(double), NULL, &ret);
	cl_v[0] = clCreateBuffer(context, CL_MEM_READ_WRITE, N * sizeof(double), NULL, &ret);
	cl_u[1] = clCreateBuffer(context, CL_MEM_READ_WRITE, N * sizeof(double), NULL, &ret);
	cl_v[1] = clCreateBuffer(context, CL_MEM_READ_WRITE, N * sizeof(double), NULL, &ret);
	cl_uhat[0] = clCreateBuffer(context, CL_MEM_READ_WRITE, N * sizeof(double), NULL, &ret);
	cl_vhat[0] = clCreateBuffer(context, CL_MEM_READ_WRITE, N * sizeof(double), NULL, &ret);
	cl_uhat[1] = clCreateBuffer(context, CL_MEM_READ_WRITE, N * sizeof(double), NULL, &ret);
	cl_vhat[1] = clCreateBuffer(context, CL_MEM_READ_WRITE, N * sizeof(double), NULL, &ret);
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
//first x
	cl_x = clCreateBuffer(context, CL_MEM_READ_WRITE, Nx * sizeof(double), NULL, &ret);
        ret = clSetKernelArg(grid, 0, sizeof(cl_mem), (void *)&cl_x);
	ret = clSetKernelArg(grid, 1, sizeof(double),(void*)&Lx);
	ret = clSetKernelArg(grid, 2, sizeof(int),(void*)&Nx);
	size_t global_work_size_x[3] = {Nx, 0, 0};
        ret = clEnqueueNDRangeKernel(command_queue, grid, 1, NULL, global_work_size_x, NULL, 0, NULL, NULL);
	ret = clFinish(command_queue);
        ret = clEnqueueReadBuffer(command_queue, cl_x, CL_TRUE, 0, Nx * sizeof(double), x, 0, NULL, NULL);
	ret = clFinish(command_queue);
//then y
	cl_y = clCreateBuffer(context, CL_MEM_READ_WRITE, Ny * sizeof(double), NULL, &ret);	
	ret = clSetKernelArg(grid, 0, sizeof(cl_mem), (void *)&cl_y);
	ret = clSetKernelArg(grid, 1, sizeof(double),(void*)&Ly);
	ret = clSetKernelArg(grid, 2, sizeof(int),(void*)&Ny);
	size_t global_work_size_y[3] = {Ny, 0, 0};

	ret = clEnqueueNDRangeKernel(command_queue, grid, 1, NULL, global_work_size_y, NULL, 0, NULL, NULL);
	ret = clFinish(command_queue);
        ret = clEnqueueReadBuffer(command_queue, cl_y, CL_TRUE, 0, Ny * sizeof(double), y, 0, NULL, NULL);
	ret = clFinish(command_queue);

//last z
	cl_z = clCreateBuffer(context, CL_MEM_READ_WRITE, Nz * sizeof(double), NULL, &ret);
	ret = clSetKernelArg(grid, 0, sizeof(cl_mem), (void *)&cl_z);
	ret = clSetKernelArg(grid, 1, sizeof(double),(void*)&Lz);
	ret = clSetKernelArg(grid, 2, sizeof(int),(void*)&Nz);
	size_t global_work_size_z[3] = {Nz, 0, 0};
	ret = clEnqueueNDRangeKernel(command_queue, grid, 1, NULL, global_work_size_z, NULL, 0, NULL, NULL);
	ret = clFinish(command_queue);
	ret = clEnqueueReadBuffer(command_queue, cl_z, CL_TRUE, 0, Nz * sizeof(double), z, 0, NULL, NULL);
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


        ret = clSetKernelArg(initialdata, 0, sizeof(cl_mem),(void *)&cl_u[0]);
	ret = clSetKernelArg(initialdata, 1, sizeof(cl_mem),(void* )&cl_v[0]);
        ret = clSetKernelArg(initialdata, 2, sizeof(cl_mem),(void *)&cl_u[1]);
	ret = clSetKernelArg(initialdata, 3, sizeof(cl_mem),(void* )&cl_v[1]);
	ret = clSetKernelArg(initialdata, 4, sizeof(cl_mem),(void* )&cl_x);
	ret = clSetKernelArg(initialdata, 5, sizeof(cl_mem),(void* )&cl_y);
	ret = clSetKernelArg(initialdata, 6, sizeof(cl_mem),(void* )&cl_z);
	ret = clSetKernelArg(initialdata, 7, sizeof(int),(void* )&Nx);
	ret = clSetKernelArg(initialdata, 8, sizeof(int),(void* )&Ny);
	ret = clSetKernelArg(initialdata, 9, sizeof(int),(void* )&Nz);
	size_t global_work_size[3] = {N, 0, 0};
        ret = clEnqueueNDRangeKernel(command_queue, initialdata, 1, NULL, global_work_size, NULL, 0, NULL, NULL);
	ret = clFinish(command_queue);
	ret = clReleaseKernel(initialdata); ret = clReleaseProgram(p_initialdata);
        ret = clEnqueueReadBuffer(command_queue, cl_u[0], CL_TRUE, 0, N * sizeof(double), u[0], 0, NULL, NULL);
	ret = clFinish(command_queue);
        ret = clEnqueueReadBuffer(command_queue, cl_v[0], CL_TRUE, 0, N * sizeof(double), v[0], 0, NULL, NULL);
	ret = clFinish(command_queue);
	ret = clReleaseMemObject(cl_x);
	ret = clReleaseMemObject(cl_y);
	ret = clReleaseMemObject(cl_z);
//write to disk
	fp=fopen("./data/xcoord.dat","w");
    	if (!fp) {fprintf(stderr, "Failed to write xcoord.dat.\n"); exit(1); }
	for(i=0;i<Nx;i++){fprintf(fp,"%lf\n",x[i]);}
    	fclose( fp );
	fp=fopen("./data/ycoord.dat","w");
    	if (!fp) {fprintf(stderr, "Failed to write ycoord.dat.\n"); exit(1); }
	for(i=0;i<Ny;i++){fprintf(fp,"%lf\n",y[i]);}
    	fclose( fp );
	fp=fopen("./data/zcoord.dat","w");
    	if (!fp) {fprintf(stderr, "Failed to write zcoord.dat.\n"); exit(1); }
	for(i=0;i<Nz;i++){fprintf(fp,"%lf\n",z[i]);}
    	fclose( fp );
	free(x); free(y); free(z);
	n=0;
	plotnum=0;
//output of initial data U
	char tmp_str[10];
	strcpy(nameconfig,"./data/u");
	sprintf(tmp_str,"%d",10000000+plotnum);
	strcat(nameconfig,tmp_str);
	strcat(nameconfig,".datbin");
	fp=fopen(nameconfig,"wb");
    	if (!fp) {fprintf(stderr, "Failed to write initialdata.\n"); exit(1); }
	for(i=0;i<N;i++){fwrite(&u[0][i], sizeof(double), 1, fp);}
    	fclose( fp );	
//V
	strcpy(nameconfig,"./data/v");
	sprintf(tmp_str,"%d",10000000+plotnum);
	strcat(nameconfig,tmp_str);
	strcat(nameconfig,".datbin");
	fp=fopen(nameconfig,"wb");
    	if (!fp) {fprintf(stderr, "Failed to write initialdata.\n"); exit(1); }
	for(i=0;i<N;i++){fwrite(&v[0][i], sizeof(double), 1, fp);}
    	fclose( fp );


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
        ret = clEnqueueNDRangeKernel(command_queue, frequencies, 1, NULL, global_work_size_x, NULL, 0, NULL, NULL);
	ret = clFinish(command_queue);
//then y
	cl_ky = clCreateBuffer(context, CL_MEM_READ_WRITE, Ny * sizeof(double), NULL, &ret);	
	ret = clSetKernelArg(frequencies, 0, sizeof(cl_mem), (void *)&cl_ky);
	ret = clSetKernelArg(frequencies, 1, sizeof(double),(void*)&Ly);
	ret = clSetKernelArg(frequencies, 2, sizeof(int),(void*)&Ny);
	ret = clEnqueueNDRangeKernel(command_queue, frequencies, 1, NULL, global_work_size_y, NULL, 0, NULL, NULL);
	ret = clFinish(command_queue);
//last z
	cl_kz = clCreateBuffer(context, CL_MEM_READ_WRITE, Nz * sizeof(double), NULL, &ret);
	ret = clSetKernelArg(frequencies, 0, sizeof(cl_mem), (void *)&cl_kz);
	ret = clSetKernelArg(frequencies, 1, sizeof(double),(void*)&Lz);
	ret = clSetKernelArg(frequencies, 2, sizeof(int),(void*)&Nz);
	ret = clEnqueueNDRangeKernel(command_queue, frequencies, 1, NULL, global_work_size_z, NULL, 0, NULL, NULL);
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

	printf("Got initial data, starting timestepping\n");
  gettimeofday(&tvs, NULL); 
	for(n=0;n<=Tmax;n++){
//linear
	ret = clfftEnqueueTransform(planHandle, CLFFT_FORWARD, 1, &command_queue, 0, NULL, NULL,cl_u, cl_uhat, tmpBufferu);
	ret = clfftEnqueueTransform(planHandle, CLFFT_FORWARD, 1, &command_queue, 0, NULL, NULL,cl_v, cl_vhat, tmpBufferv);
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
	ret = clSetKernelArg(linearpart, 12, sizeof(double),(void* )&b[0]);
	ret = clSetKernelArg(linearpart, 13, sizeof(double),(void* )&b[1]);
	ret = clSetKernelArg(linearpart, 14, sizeof(int),(void* )&Nx);
	ret = clSetKernelArg(linearpart, 15, sizeof(int),(void* )&Ny);
	ret = clSetKernelArg(linearpart, 16, sizeof(int),(void* )&Nz);
        ret = clEnqueueNDRangeKernel(command_queue, linearpart, 1, NULL, global_work_size, NULL, 0, NULL, NULL);
	ret = clFinish(command_queue);

	ret = clfftEnqueueTransform(planHandle, CLFFT_BACKWARD, 1, &command_queue, 0, NULL, NULL,cl_uhat, cl_u, tmpBufferu);
	ret = clfftEnqueueTransform(planHandle, CLFFT_BACKWARD, 1, &command_queue, 0, NULL, NULL,cl_vhat, cl_v, tmpBufferv);
	ret = clFinish(command_queue);    
//nonlinearpart
        ret = clSetKernelArg(nonlinearpart, 0, sizeof(cl_mem),(void *)&cl_u[0]);
        ret = clSetKernelArg(nonlinearpart, 1, sizeof(cl_mem),(void *)&cl_u[1]);
	ret = clSetKernelArg(nonlinearpart, 2, sizeof(cl_mem),(void* )&cl_v[0]);
	ret = clSetKernelArg(nonlinearpart, 3, sizeof(cl_mem),(void* )&cl_v[1]);
	ret = clSetKernelArg(nonlinearpart, 4, sizeof(double),(void* )&dt);
	ret = clSetKernelArg(nonlinearpart, 5, sizeof(double),(void* )&a[0]);
	ret = clSetKernelArg(nonlinearpart, 6, sizeof(double),(void* )&a[1]);
        ret = clEnqueueNDRangeKernel(command_queue, nonlinearpart, 1, NULL, global_work_size, NULL, 0, NULL, NULL);
	ret = clFinish(command_queue);		
// linear part
	ret = clfftEnqueueTransform(planHandle, CLFFT_FORWARD, 1, &command_queue, 0, NULL, NULL,cl_u, cl_uhat, tmpBufferu);
	ret = clfftEnqueueTransform(planHandle, CLFFT_FORWARD, 1, &command_queue, 0, NULL, NULL,cl_v, cl_vhat, tmpBufferv);	
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
	ret = clSetKernelArg(linearpart, 12, sizeof(double),(void* )&b[0]);
	ret = clSetKernelArg(linearpart, 13, sizeof(double),(void* )&b[1]);
	ret = clSetKernelArg(linearpart, 14, sizeof(int),(void* )&Nx);
	ret = clSetKernelArg(linearpart, 15, sizeof(int),(void* )&Ny);
	ret = clSetKernelArg(linearpart, 16, sizeof(int),(void* )&Nz);
        ret = clEnqueueNDRangeKernel(command_queue, linearpart, 1, NULL, global_work_size, NULL, 0, NULL, NULL);
	ret = clFinish(command_queue);

	ret = clfftEnqueueTransform(planHandle, CLFFT_BACKWARD, 1, &command_queue, 0, NULL, NULL,cl_uhat, cl_u, tmpBufferu);
	ret = clfftEnqueueTransform(planHandle, CLFFT_BACKWARD, 1, &command_queue, 0, NULL, NULL,cl_vhat, cl_v, tmpBufferv);
	ret = clFinish(command_queue);
// done
	if(n==plottime){
		printf("time:%lf, step:%d,%d\n",n*dt,n,plotnum);
		plottime=plottime+plotgap;
		plotnum=plotnum+1;
        ret = clEnqueueReadBuffer(command_queue, cl_u[0], CL_TRUE, 0, N * sizeof(double), u[0], 0, NULL, NULL);
        ret = clEnqueueReadBuffer(command_queue, cl_v[0], CL_TRUE, 0, N * sizeof(double), v[0], 0, NULL, NULL);
	ret = clFinish(command_queue);
//output of data U
	char tmp_str[10];
	strcpy(nameconfig,"./data/u");
	sprintf(tmp_str,"%d",10000000+plotnum);
	strcat(nameconfig,tmp_str);
	strcat(nameconfig,".datbin");
	fp=fopen(nameconfig,"wb");
    	if (!fp) {fprintf(stderr, "Failed to write u-data.\n"); exit(1); }
	for(i=0;i<N;i++){fwrite(&u[0][i], sizeof(double), 1, fp);}
    	fclose( fp );	
//V
	strcpy(nameconfig,"./data/v");
	sprintf(tmp_str,"%d",10000000+plotnum);
	strcat(nameconfig,tmp_str);
	strcat(nameconfig,".datbin");
	fp=fopen(nameconfig,"wb");
    	if (!fp) {fprintf(stderr, "Failed to write v-data.\n"); exit(1); }
	for(i=0;i<N;i++){fwrite(&v[0][i], sizeof(double), 1, fp);}
    	fclose( fp );
}
	}
 	gettimeofday(&tve, NULL); 
	printf("Finished time stepping\n");
 	elapsedTime = (tve.tv_sec - tvs.tv_sec) * 1000.0;      // sec to ms
    	elapsedTime += (tve.tv_usec - tvs.tv_usec) / 1000.0;   // us to ms
   	printf("Program took %lfms for execution\n",elapsedTime);



	clReleaseMemObject(cl_u[0]);
	clReleaseMemObject(cl_u[1]);
	clReleaseMemObject(cl_v[0]);
	clReleaseMemObject(cl_v[1]);
	clReleaseMemObject(cl_uhat[0]);
	clReleaseMemObject(cl_uhat[1]);
	clReleaseMemObject(cl_vhat[0]);
	clReleaseMemObject(cl_vhat[1]);
	clReleaseMemObject(cl_kx);
	clReleaseMemObject(cl_ky);
	clReleaseMemObject(cl_kz);
    	ret = clReleaseKernel(frequencies); ret = clReleaseProgram(p_frequencies);
    	ret = clReleaseKernel(linearpart); ret = clReleaseProgram(p_linearpart);
    	ret = clReleaseKernel(nonlinearpart); ret = clReleaseProgram(p_nonlinearpart);
	free(u[0]);
	free(v[0]);
	clReleaseMemObject(tmpBufferu);
	clReleaseMemObject(tmpBufferv);
	/* Release the plan. */
	ret = clfftDestroyPlan(&planHandle);
	/* Release clFFT library. */
	clfftTeardown();

	ret = clReleaseCommandQueue(command_queue);
     	ret = clReleaseContext(context);	
	printf("Program execution complete\n");

	return 0;
}
