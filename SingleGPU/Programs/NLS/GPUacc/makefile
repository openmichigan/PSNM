F90   = pgf90
FLAGS = -O3 -acc -Mcuda -ta=nvidia -Minfo=accel
CUDA_PATH = /usr/local/pgi/linux86-64/2012/cuda/4.1
LIBS = -I$(CUDA_PATH)/include -L$(CUDA_PATH)/lib64 -lcudart -lcufft

BIN =  cubicNLS_acc

all: $(BIN)

cubicNLS_acc: cubicNLS.f90 
	$(F90) $(FLAGS) -o $@ $^ $(LIBS)

clean:
	$(RM) $(BIN) *.o *.mod
     
