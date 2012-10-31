NVCC=nvcc
NVCFLAGS=-m64 -O3 -arch=sm_20
NVLIBS=-L/uf/ac/rigge/NVIDIA_GPU_Computing_SDK/CUDALibraries/common/lib/ -lcutil_x86_64 -lcufft $(FFTW33_LIBPATH) -lfftw3
NVINC=-I/uf/ac/rigge/NVIDIA_GPU_Computing_SDK/CUDALibraries/common/inc/ $(FFTW33_INCPATH)
F90=pgf90
CUDA_PATH=/usr/local/pgi/linux86-64/2012/cuda/4.1/
CCFLAGS=-I$(CUDA_PATH)/include -L$(CUDA_PATH)/lib64 -lcudart -lcufft
SRCS=enercalc.f90 getgrid.f90 initialdata.f90 savedata.f90 saveresults.f90 sgsemiimp2d.f90 storeold.f90 

all: sg2d
#sg2df

enercalc.o: enercalc.f90
	gfortran -static-libgfortran -O3 -m64 -c enercalc.f90

sg2d: sg2d.cu enercalc.o
	$(NVCC) $(NVCFLAGS) $(NVLIBS) $(NVINC) sg2d.cu -c
	$(NVCC) $(NVCFLAGS) $(NVLIBS) sg2d.o enercalc.o -o sg2d -lgfortran -g

sg2df: $(SRCS)
	$(F90) $(OMPFLAGS) $(CCFLAGS) -o $@ $(SRCS)

clean:
	rm sg2d *.o
#*.mod
