F90       = pgf90
FLAGS = -O3 -acc -Mcuda -ta=nvidia
CUDA_PATH = /usr/local/pgi/linux86-64/2012/cuda/4.1
CCFLAGS  = -I$(CUDA_PATH)/include -L$(CUDA_PATH)/lib64 -lcudart -lcufft

BIN =  ns2dcn_acc

all: $(BIN)

ns2dcn_acc: navierstokes.f90 
	$(F90) $(FLAGS) $(CCFLAGS) -o $@ $<

clean:
	$(RM) *.mod
	$(RM) $(BIN)
