F90       = pgf90
FLAGS = -O3 -Mcuda -Minfo 
CUDA_PATH = /usr/local/pgi/linux86-64/2012/cuda/4.1
CCFLAGS  = -I$(CUDA_PATH)/include -L$(CUDA_PATH)/lib64 -lcudart -lcufft

BIN =  ns2dcn

all: $(BIN)

ns2dcn: navierstokes.cuf 
	$(F90) $(FLAGS) $(CCFLAGS) -o $@ $<

clean:
	$(RM) *.mod
	$(RM) $(BIN)
