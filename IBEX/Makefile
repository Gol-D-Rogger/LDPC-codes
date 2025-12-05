#Makefile for IBEX
CC = g++
CXX = g++

# Target executable
TARGET_EXEC ?= ssd_fc

BUILD_DIR ?= ./build
SRC_DIRS ?= ./src
DUMPS_DIR ?= ./output

# Common source files
COMMON_C_SRCS := $(shell find $(SRC_DIRS) -name '*.c' ! -name '*.bak')
COMMON_CPP_SRCS := $(SRC_DIRS)/fc_dsp.cpp $(SRC_DIRS)/transceiver.cpp $(SRC_DIRS)/isoc23_compat.cpp

# Main program sources
MAIN_SRCS := $(SRC_DIRS)/SSD_FC.cpp $(SRC_DIRS)/ldpc_codec.cpp

# Object files
COMMON_C_OBJS := $(COMMON_C_SRCS:$(SRC_DIRS)/%.c=$(BUILD_DIR)/$(SRC_DIRS)/%.c.o)
COMMON_CPP_OBJS := $(COMMON_CPP_SRCS:$(SRC_DIRS)/%.cpp=$(BUILD_DIR)/$(SRC_DIRS)/%.cpp.o)
MAIN_OBJS := $(MAIN_SRCS:$(SRC_DIRS)/%.cpp=$(BUILD_DIR)/$(SRC_DIRS)/%.cpp.o)
OBJS := $(COMMON_C_OBJS) $(COMMON_CPP_OBJS) $(MAIN_OBJS)

INC_DIRS := $(shell find $(SRC_DIRS) -type d)
INC_FLAGS := $(addprefix -I,$(INC_DIRS))

CPPFLAGS ?= $(INC_FLAGS) -MMD -MP -g
LDFLAGS ?= -lm

# Default target: build the program
.PHONY: all
all: $(TARGET_EXEC)

# Build executable
$(TARGET_EXEC): $(OBJS)
	$(CXX) $(OBJS) -o $@ $(LDFLAGS)

# C source compilation
$(BUILD_DIR)/$(SRC_DIRS)/%.c.o: $(SRC_DIRS)/%.c
	$(MKDIR_P) $(dir $@)
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -g -o $@

# C++ source compilation
$(BUILD_DIR)/$(SRC_DIRS)/%.cpp.o: $(SRC_DIRS)/%.cpp
	$(MKDIR_P) $(dir $@)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -g -o $@

.PHONY: clean
clean:
	$(RM) -r $(BUILD_DIR) $(TARGET_EXEC)
	$(RM) $(DUMPS_DIR)/*_matrix.txt
	$(RM) $(DUMPS_DIR)/*_matrix_*.txt

.PHONY: run
run:
	./$(TARGET_EXEC) LDPC config/Ibex_hd_row8.cnfg AWGN 5.30 | tee perf/test.log

MKDIR_P ?= mkdir -p

-include $(OBJS:.o=.d)
