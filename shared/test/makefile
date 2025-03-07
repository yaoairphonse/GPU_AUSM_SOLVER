CUDA_LIB := /usr/local/cuda/lib64 -lcuda -lcudart

all: CPU GPU
	g++ main.o cuda.o -o song.run -L $(CUDA_LIB)

CPU:
	g++ main.cpp -c

GPU:
	nvcc cuda.cu -c -arch sm_35

