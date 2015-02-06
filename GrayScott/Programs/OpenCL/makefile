####### Compiler, tools and options

CC            = gcc
CFLAGS        = -m64 -pipe -O3 -Wno-unused-but-set-parameter -Wall
DEL_FILE      = rm -f

LIBS = -lOpenCL -I/home/quell/grayscott/clFFT/src/package/include/ -L/home/quell/grayscott/clFFT/src/package/lib64 -lclFFT -lm
 
#

####### Build rules 

CL: main_gs.c
	$(DEL_FILE) gs
	$(DEL_FILE) *.o
	$(CC) $(CFLAGS) -o gs main_gs.c $(LIBS)
	export LD_LIBRARY_PATH=/home/quell/grayscott/clFFT/src/package/lib64:$$LD_LIBRARY_PATH




out: xdmfcreate.f90
	gfortran -o out xdmfcreate.f90


clean: 

	$(DEL_FILE) gs
	$(DEL_FILE) out
	$(DEL_FILE) *.xmf
	$(DEL_FILE) *.o
	$(DEL_FILE) ./data/u*
	$(DEL_FILE) ./data/v*


cleandata: 

	$(DEL_FILE) ./data/u*
	$(DEL_FILE) ./data/v*
