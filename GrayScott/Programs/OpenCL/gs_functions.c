//
//
//
//This file contains only functions for main_gs.c
//
//
#include <CL/cl.h>
#include <sys/time.h>
#include <string.h>

//read the INPUTFILE
void parainit(int * Nx, int * Ny, int * Tmax, int * plotgap, double * Lx, double * Ly, double * dt, double * Du, double * Dv, double * A, double *B ){

	int intcomm[4];
	double dpcomm[7];
		char	InputFileName[]="./INPUTFILE";
		FILE*fp;
		fp=fopen(InputFileName,"r");
   		 if(!fp) {fprintf(stderr, "Failed to load IPUTFILE.\n");exit(1);}	 
		int ierr=fscanf(fp, "%d %d %d %d %lf %lf %lf %lf %lf %lf %lf", &intcomm[0],&intcomm[1],&intcomm[2],&intcomm[3],&dpcomm[0],&dpcomm[1],&dpcomm[2],&dpcomm[3],&dpcomm[4],&dpcomm[5],&dpcomm[6]);
		if(ierr!=11){fprintf(stderr, "INPUTFILE corrupted:%d\n",ierr);exit(1);}	
		fclose(fp);
		printf("NX %d\nNY %d\nTmax %d\nplotgap %d\n",intcomm[0],intcomm[1],intcomm[2],intcomm[3]); 
		printf("Lx %lf\nLy %lf\ndt %lf\nDu %lf\nDv %lf\nF %lf\nk %lf\n",dpcomm[0],dpcomm[1],dpcomm[2],dpcomm[3],dpcomm[4],dpcomm[5],dpcomm[6]);
	*Nx=intcomm[0];
	*Ny=intcomm[1];
	*Tmax=intcomm[2];
	*plotgap=intcomm[3];
	*Lx=dpcomm[0];
	*Ly=dpcomm[1];
	*dt=dpcomm[2];
	*Du=dpcomm[3];
	*Dv=dpcomm[4];
	*A=dpcomm[5];
	*B=dpcomm[5]+dpcomm[6];
	printf("Read Inputfile\n");
};


//loads a kernel from a file
#define MAX_SOURCE_SIZE 8192
void loadKernel(cl_kernel *kernel,cl_context *context, cl_device_id *device_id, char*name){
        cl_program p_kernel;
        cl_int ret=0;
        size_t source_size;
        char *source_str;
        char nameconfig[100];
		int i=0;
        source_str = (char *)malloc(MAX_SOURCE_SIZE*sizeof(char));
        for(i=0;i<MAX_SOURCE_SIZE;i++){source_str[i]='\0';}
                FILE* fp;
                strcpy(nameconfig,"./");
                strcat(nameconfig,name);
                strcat(nameconfig,".cl");
                fp = fopen(nameconfig, "r");
                if (!fp) {fprintf(stderr, "Failed to load kernel.\n"); exit(1); }
                source_size = fread( source_str, sizeof(char), MAX_SOURCE_SIZE, fp );
                fclose( fp );

        p_kernel = clCreateProgramWithSource(*context, 1, (const char **)&source_str, (const size_t *)&source_size, &ret);
	if(ret!=CL_SUCCESS){printf("createProgram ret:%d\n",ret);exit(1); }
        ret = clBuildProgram(p_kernel, 1, &*device_id, NULL, NULL, NULL);
        if(ret!=CL_SUCCESS){printf("buildProgram ret:%d\n",ret); exit(1); }
        *kernel = clCreateKernel(p_kernel, name, &ret);
        if(ret!=CL_SUCCESS){printf("createKernel ret:%d\n",ret);exit(1); }
        ret = clReleaseProgram(p_kernel);
        if(ret!=CL_SUCCESS){printf("releaseProgram ret:%d\n",ret);exit(1); }
		printf("got kernel %s\n",name);
        free(source_str);
};

//displays an array on gpu memory (debug)
void printCL(cl_mem* cl_u,cl_command_queue* command_queue, int Nx,int Ny){
	double* u;
	int i=0;
	int j=0;
	cl_int ret=0;
	u=(double*)malloc(Nx*Ny*sizeof(double));
  	ret = clEnqueueReadBuffer(*command_queue, *cl_u, CL_TRUE, 0, Nx*Ny*sizeof(double), u, 0, NULL, NULL);
	ret = clFinish(*command_queue); if(ret!=CL_SUCCESS){printf("failed");}
	for(i=0;i<Nx;i++){
		for(j=0;j<Ny;j++){
			printf("%f ",u[i+Nx*j]);
		}
		printf("\n");
	}	
	printf("\n");
	free(u);
};

//displays an real part of complex array on gpu memory (debug)
void printCL_C(cl_mem* cl_u,cl_command_queue* command_queue, int Nx,int Ny){
	double* u;
	int i=0;
	int j=0;
	cl_int ret=0;
	u=(double*)malloc(2*Nx*Ny*sizeof(double));
  	ret = clEnqueueReadBuffer(*command_queue, *cl_u, CL_TRUE, 0, 2*Nx*Ny*sizeof(double), u, 0, NULL, NULL);
	ret = clFinish(*command_queue);
	if(ret!=CL_SUCCESS){printf("failed");}
	
	for(i=0;i<Nx;i++){
		for(j=0;j<Ny;j++){
			printf("%f ",u[2*i+Nx*2*j]);
		}
		printf("\n");
	}	
	printf("\n");
	free(u);
};

//make plans for FFT
void fftinit(clfftPlanHandle *planHandle, cl_context* context, cl_command_queue* command_queue,	cl_mem* tmpBuffer, int Nx,int Ny){
	clfftDim dim = CLFFT_2D;
	size_t clLength[2] = {Nx,Ny};
	cl_int ret=0;

	// Setup clFFT. 
	clfftSetupData fftSetup;
	ret = clfftInitSetupData(&fftSetup);
	if(ret!=CL_SUCCESS){printf("clFFT init ret:%d\n",ret);exit(1); }
	ret = clfftSetup(&fftSetup);
	if(ret!=CL_SUCCESS){printf("clFFT Setup ret:%d\n",ret);exit(1); }
	// Create a default plan for a complex FFT. 
	ret = clfftCreateDefaultPlan(&*planHandle, *context, dim, clLength);
	if(ret!=CL_SUCCESS){printf("clFFT Plan ret:%d\n",ret);exit(1); }
	// Set plan parameters. 
	ret = clfftSetPlanPrecision(*planHandle, CLFFT_DOUBLE);
	if(ret!=CL_SUCCESS){printf("clFFT Precision ret:%d\n",ret);exit(1); }
	//ret = clfftSetPlanBatchSize(*planHandle, (size_t) Ny );
	//if(ret!=CL_SUCCESS){printf("clFFT Batch ret:%d\n",ret);exit(1); }
	ret = clfftSetLayout(*planHandle, CLFFT_COMPLEX_INTERLEAVED, CLFFT_COMPLEX_INTERLEAVED);
	if(ret!=CL_SUCCESS){printf("clFFT Layout ret:%d\n",ret);exit(1); }
	ret = clfftSetResultLocation(*planHandle, CLFFT_OUTOFPLACE);
	if(ret!=CL_SUCCESS){printf("clFFT Place ret:%d\n",ret);exit(1); }

	// Bake the plan. 
	ret = clfftBakePlan(*planHandle, 1, &*command_queue, NULL, NULL);
	if(ret!=CL_SUCCESS){printf("clFFT Bake ret:%d\n",ret);exit(1); }
	// Create temporary buffer. 
	// Size of temp buffer. 
	size_t tmpBufferSize = 0;
	ret = clfftGetTmpBufSize(*planHandle, &tmpBufferSize);
	if ((ret == CL_SUCCESS) && (tmpBufferSize > 0)) {
		*tmpBuffer = clCreateBuffer(*context, CL_MEM_READ_WRITE, tmpBufferSize, NULL, &ret);
		if (ret != CL_SUCCESS){printf("Error with tmpBuffer clCreateBuffer\n");exit(1);}
	}

};

//destroy plans
void fftdestroy(clfftPlanHandle *planHandle,cl_mem* tmpBuffer){
	cl_int ret=0;	
	clReleaseMemObject(*tmpBuffer);
	ret = clfftDestroyPlan(&*planHandle);
	if(ret!=0){printf("Error while destroying fft");exit(1);}
	clfftTeardown();
};

//fft2dfoward
void fft2dfor(cl_mem *cl_u, cl_mem *cl_uhat, clfftPlanHandle *planHandle, cl_command_queue* command_queue, cl_mem* tmpBuffer){
	int ret=0;
	ret = clfftEnqueueTransform(*planHandle, CLFFT_FORWARD, 1, command_queue, 0, NULL, NULL,&*cl_u, &*cl_uhat, *tmpBuffer);
	if (ret != CL_SUCCESS){printf("FFT failedA%d",ret);}
	ret = clFinish(*command_queue);
	if (ret != CL_SUCCESS){printf("FFT failedB%d",ret);}
};

//fft2dback
void fft2dback(cl_mem *cl_u, cl_mem *cl_uhat, clfftPlanHandle *planHandle, cl_command_queue* command_queue, cl_mem* tmpBuffer){
	int ret=0;
	ret = clfftEnqueueTransform(*planHandle, CLFFT_BACKWARD, 1, command_queue, 0, NULL, NULL,&*cl_uhat, &*cl_u, *tmpBuffer);
	if (ret != CL_SUCCESS){printf("FFT failedC%d",ret);}
	ret = clFinish(*command_queue);	
	if (ret != CL_SUCCESS){printf("FFT failedD%d",ret);}
};

		
//writes an image to disk and returns the maximum of cl_u
double writeimage(cl_mem* cl_u, cl_command_queue *command_queue, int Nx,int Ny, int plotnum, char* prefix){
	int i=0;
    cl_int ret=0;
    int header=54;
	double* u;
    u=(double*)malloc(2*Nx*Ny*sizeof(double));
    ret = clEnqueueReadBuffer(*command_queue, *cl_u, CL_TRUE, 0, 2*Nx*Ny * sizeof(double), u, 0, NULL, NULL);
    ret = clFinish(*command_queue);
    if(ret!=0){printf("Error hahah");}
    double max=0.0;
    for(i=0;i<Nx*Ny;i++){
     	if(u[2*i]>max){max=u[2*i];}
    }
    unsigned char*picture=(unsigned char*)malloc((3*Nx*Ny+header)*sizeof(unsigned char));

    for(i=0;i<Nx*Ny;i++){
        picture[3*i+header+0]=(unsigned char)(255*u[2*i]/max);
        picture[3*i+header+1]=(unsigned char)(255*u[2*i]/max);
        picture[3*i+header+2]=(unsigned char)(255*u[2*i]/max);
    }
	//header for bmp file
	int w=Ny;
	int h=Nx;
	int padSize=(4-w%4)%4;
	int filesize=header + 3*h*w+h*padSize;
	unsigned char bmppad[3] = {0,0,0}; //padding
	picture[ 0]='B';
	picture[ 1]='M';
	picture[ 2] = (unsigned char)(filesize    );
	picture[ 3] = (unsigned char)(filesize>> 8);
	picture[ 4] = (unsigned char)(filesize>>16);
	picture[ 5] = (unsigned char)(filesize>>24);
	picture[ 6] = 0;
	picture[ 7] = 0;
	picture[ 8] = 0;
	picture[ 9] = 0;
	picture[10] = 54;
	picture[11] = 0;
	picture[12] = 0;
	picture[13] = 0;
	picture[14] = 40;
	picture[15] = 0;
	picture[16] = 0;
	picture[17] = 0;//3
	picture[18] = (unsigned char)(       w    );
	picture[19] = (unsigned char)(       w>> 8);
	picture[20] = (unsigned char)(       w>>16);
	picture[21] = (unsigned char)(       w>>24);
	picture[22] = (unsigned char)(       h    );
	picture[23] = (unsigned char)(       h>> 8);
	picture[24] = (unsigned char)(       h>>16);
	picture[25] = (unsigned char)(       h>>24);
	picture[26] = 1;
	picture[27] = 0;
	picture[28] = 24;
	picture[29] = 0;
	for(i=30;i<54;i++){
		picture[i]=0;
	}
	FILE*fp;
	//file name
	char tmp_str[10];
	char nameconfig[100];
	strcpy(nameconfig,"./data/");
	strcat(nameconfig,prefix);
	sprintf(tmp_str,"%d",10000000+plotnum);
	strcat(nameconfig,tmp_str);
	strcat(nameconfig,".bmp");
	fp=fopen(nameconfig,"wb");
    if (!fp) {fprintf(stderr, "Failed to write data.\n"); exit(1); }
	for(i=0;i<header;i++){fwrite(&picture[i], sizeof(unsigned char), 1, fp);}
	for(i=0;i<h;i++){
		fwrite(picture+(w*(h-i-1)*3)+header,3* sizeof(unsigned char), w, fp);
		fwrite(bmppad,sizeof(unsigned char),(4-(w*3)%4)%4,fp);
	}
    fclose( fp );
	free(picture);
	free(u);
	return max;	
};

//writes the array to disk (debug)
void writedata(cl_mem* cl_u, cl_command_queue *command_queue, int Nx,int Ny, int plotnum,char* prefix){
	int i=0;
    cl_int ret=0;
	double* u;
    u=(double*)malloc(Nx*Ny*sizeof(double));
    ret = clEnqueueReadBuffer(*command_queue, *cl_u, CL_TRUE, 0, Nx*Ny * sizeof(double), u, 0, NULL, NULL);
    ret = clFinish(*command_queue);
    if(ret!=0){printf("Error hahah");}

	FILE*fp;
	//file name
	char tmp_str[10];
	char nameconfig[100];
	strcpy(nameconfig,"./data/");
	strcat(nameconfig,prefix);
	sprintf(tmp_str,"%d",10000000+plotnum);
	strcat(nameconfig,tmp_str);
	strcat(nameconfig,".datbin");
	fp=fopen(nameconfig,"wb");
    if (!fp) {fprintf(stderr, "Failed to write data.\n"); exit(1); }
	for(i=0;i<Nx;i++){fwrite(u+i*Ny, sizeof(double), Ny, fp);}
    fclose( fp );
	free(u);
};

//writes the real part of complex array to disk
void writedata_C(cl_mem* cl_u, cl_command_queue *command_queue, int Nx,int Ny, int plotnum,char* prefix){
	int i=0;
    cl_int ret=0;
	double* u;
    u=(double*)malloc(2*Nx*Ny*sizeof(double));
    ret = clEnqueueReadBuffer(*command_queue, *cl_u, CL_TRUE, 0, 2*Nx*Ny * sizeof(double), u, 0, NULL, NULL);
    ret = clFinish(*command_queue);
    if(ret!=0){printf("Error hahah");}
	for(i=0;i<Nx*Ny;i++){u[i]=u[2*i];}
	FILE*fp;
	//file name
	char tmp_str[10];
	char nameconfig[100];
	strcpy(nameconfig,"./data/");
	strcat(nameconfig,prefix);
	sprintf(tmp_str,"%d",10000000+plotnum);
	strcat(nameconfig,tmp_str);
	strcat(nameconfig,".datbin");
	fp=fopen(nameconfig,"wb");
    if (!fp) {fprintf(stderr, "Failed to write data.\n"); exit(1); }
	for(i=0;i<Nx;i++){fwrite(u+i*Ny, sizeof(double), Ny, fp);}
    fclose( fp );
	free(u);
};

//loades the data from disk (debug)
void loaddata(cl_mem* cl_u, cl_command_queue *command_queue, int Nx,int Ny, char*name){
	int i=0;
	int numread=0;
    cl_int ret=0;
	double* u;
    u=(double*)malloc(Nx*Ny*sizeof(double));
	FILE*fp;
	fp=fopen(name,"rb");
    if (!fp) {fprintf(stderr, "Failed to open file.\n"); exit(1); }
	for(i=0;i<Nx;i++){
		numread=fread(u+i*Ny, sizeof(double), Ny, fp);
    if (numread!=Ny) {fprintf(stderr, "Failed to read file.\n"); exit(1); }
	}
    fclose( fp );

    ret = clEnqueueWriteBuffer(*command_queue, *cl_u, CL_TRUE, 0, Nx*Ny * sizeof(double), u, 0, NULL, NULL);
    ret = clFinish(*command_queue);
    if(ret!=0){printf("Error hahah");}

	free(u);
};

//writes an array to disc ASCII 
void writearray(double*u,int lenght,char* prefix){
	int i=0;
	char nameconfig[100];
	strcpy(nameconfig,"./data/");
	strcat(nameconfig,prefix);
	FILE*fp;
	fp=fopen(nameconfig,"w");
	if (!fp) {fprintf(stderr, "Failed to write data.\n"); exit(1); }
	for(i=0;i<lenght;i++){
		fprintf(fp,"%.17g\n",u[i]);
	}
};


//start measuring time
void mtime_s(struct timeval* tvs){
 	gettimeofday(&*tvs, NULL);
};

//end measuring time+ printing
void mtime_e(struct timeval* tvs, char* printme){
	struct timeval tve;
   	double elapsedTime;
 	gettimeofday(&tve, NULL); 
 	elapsedTime = (tve.tv_sec - (*tvs).tv_sec) * 1000.0;      // sec to ms
    elapsedTime += (tve.tv_usec - (*tvs).tv_usec) / 1000.0;   // us to ms
   	printf("%s%lfms\n",printme,elapsedTime);
};
