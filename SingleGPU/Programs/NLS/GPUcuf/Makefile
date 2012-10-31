F90       = pgf90
OMPFLAGS  = -O3 -Mcuda -Minfo=accel -Mprof=lines -Mlarge_arrays
CUDA_PATH = /usr/local/pgi/linux86-64/2012/cuda/4.1
CCFLAGS   = -I$(CUDA_PATH)/include -L$(CUDA_PATH)/lib64 -lcudart -lcufft

BIN =  nls2d

all: $(BIN)

nls2d: cubicNLS.cuf 
	$(F90) $(OMPFLAGS) $(CCFLAGS) -o $@ $<

clean:
	$(RM) *.mod
	$(RM) $(BIN)
