####### Compiler, tools and options

CC            = gcc
CFLAGS        = -m64 -pipe -O2 -Wno-unused-but-set-parameter -Wall
#OPENCL        = /opt/AMDAPP
#CLFFT         = /opt/clFFT
#INC_OPENCL    = -I$(OPENCL)/include 
#INC_CLFFT     = -I$(CLFFT)/include
#INCPATH       = $(INC_OPENCL) $(INC_CLFFT) 
#LIBS_CLFFT    = -L$(CLFFT)/lib64 -lclFFT -L$(OPENCL)/lib/x86_64 -lOpenCL 
DEL_FILE      = rm -f

LIBS = -lOpenCL -I/home/quell/grayscott/clFFT/src/package/include/ -L/home/quell/grayscott/clFFT/src/package/lib64 -lclFFT -lm
 
#

####### Build rules 

CL: grayscottOpenCLs.c
	$(DEL_FILE) grayscottOpenCLs
	$(DEL_FILE) *.o
	$(CC) $(CFLAGS) -o grayscottOpenCLs grayscottOpenCLs.c $(LIBS)
	export LD_LIBRARY_PATH=/home/quell/grayscott/clFFT/src/package/lib64:$$LD_LIBRARY_PATH
	export DISPLAY=:0

CL43: grayscottOpenCL43.c
	$(DEL_FILE) grayscottOpenCL43
	$(DEL_FILE) *.o
	$(CC) $(CFLAGS) -o grayscottOpenCL43 grayscottOpenCL43.c $(LIBS)
	export LD_LIBRARY_PATH=/home/quell/grayscott/clFFT/src/package/lib64:$$LD_LIBRARY_PATH

CL4k: grayscottOpenCL4k.c
	$(DEL_FILE) grayscottOpenCL4k
	$(DEL_FILE) *.o
	$(CC) $(CFLAGS) -o grayscottOpenCL4k grayscottOpenCL4k.c $(LIBS)
	export LD_LIBRARY_PATH=/home/quell/grayscott/clFFT/src/package/lib64:$$LD_LIBRARY_PATH

info: infocl.c
	$(CC) $(CFLAGS) -o infocl infocl.c $(LIBS)

out: xdmfcreate.f90
	gfortran -o out xdmfcreate.f90

outa: xdmfcreatea.f90
	gfortran -o outa xdmfcreatea.f90
clean: 
	$(DEL_FILE) grayscottOpenCL43
	$(DEL_FILE) grayscottOpenCLs
	$(DEL_FILE) infocl
	$(DEL_FILE) out
	$(DEL_FILE) outa
	$(DEL_FILE) *.xmf
	$(DEL_FILE) *.o
