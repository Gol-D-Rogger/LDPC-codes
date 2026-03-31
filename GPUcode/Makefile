NVCC := nvcc
NVCCFLAGS := -std=c++11 -O3 -Xptxas=-v -lineinfo
GXX := g++
GXXFLAGS := -std=c++17 -O2

TARGET := a.out
CPU_TARGET := cpu.out

CPU_SRCS := \
	codeword.cpp \
	h_matrix.cpp \
	error_injector.cpp \
	codec_support.cpp \
	decoder.cpp \
	main.cpp \
	calc_llr.cpp

CUDA_SRCS := \
	$(CPU_SRCS) \
	decoder_gpu.cu

.PHONY: all cpu_build cuda_build run debug clean

all: cpu_build

cpu_build:
	$(GXX) $(GXXFLAGS) $(CPU_SRCS) -o $(CPU_TARGET)

cuda_build:
	$(NVCC) $(NVCCFLAGS) $(CUDA_SRCS) -o $(TARGET)

run: cpu_build
	./$(CPU_TARGET)

debug:
	$(GXX) -std=c++17 -O0 -g $(CPU_SRCS) -o $(CPU_TARGET)
	gdb ./$(CPU_TARGET)

clean:
	rm -f $(TARGET) $(CPU_TARGET)
