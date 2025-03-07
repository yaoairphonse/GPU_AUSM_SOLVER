CUDA_LIB := /usr/local/cuda/lib64 -lcuda -lcudart

all: CPU GPU
	g++ main.o cudaX37.o -o song.run -L $(CUDA_LIB)

CPU:
	g++ main.cpp -c

GPU:
	nvcc cudaX37.cu -c -arch sm_86

