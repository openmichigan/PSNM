CC       = nvcc
CUDA_PATH = /usr/local/cuda/
CUDA_SDK_PATH = /usr/local/NVIDIA_GPU_Computing_SDK/CUDALibraries
CCFLAGS  = -O3 -arch=sm_20 -I$(CUDA_SDK_PATH)/common/inc \
	-I$(CUDA_PATH)/include -L$(CUDA_PATH)/lib64 -lcudart -lcufft

BIN = nls2d 

all: $(BIN)

nls2d: NLSsplitting.cu
	$(CC) $(CCFLAGS) -o $@ $<
	
clean:
	$(RM) $(BIN)
