F90       = pgf90
FLAGS = -O3 -Mcuda -Minfo
CUDA_PATH = /usr/local/pgi/linux86-64/2011/cuda/4.0/
CCFLAGS  = -I$(CUDA_PATH)/include -L$(CUDA_PATH)/lib64 -lcudart -lcufft $(FFTW33_LIBPATH) -lfftw3

BIN =  sg2d

SOURCES = sgsemiimp2d.cuf getgrid.f90 initialdata.f90 enercalc.f90 savedata.f90 saveresults.f90
all: $(BIN)

sg2d: $(SOURCES)
	$(F90) -o $(BIN) $(FLAGS) $(SOURCES) $(CCFLAGS)

clean:
	$(RM) *.mod
	$(RM) *.o
	$(RM) $(BIN)
