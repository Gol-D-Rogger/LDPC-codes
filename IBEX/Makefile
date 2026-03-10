#Makefile for IBEX
CC = g++
CXX = g++

# Target executable
TARGET_EXEC ?= ssd_fc

BUILD_DIR ?= ./build
SRC_DIRS ?= ./src
SRC2_DIR ?= ./src2
DUMPS_DIR ?= ./output

# Common source files
COMMON_C_SRCS := $(shell find $(SRC_DIRS) -name '*.c' ! -name '*.bak')
COMMON_CPP_SRCS := $(SRC_DIRS)/fc_dsp.cpp $(SRC_DIRS)/transceiver.cpp $(SRC_DIRS)/isoc23_compat.cpp
SRC2_COMMON_CPP_SRCS := $(COMMON_CPP_SRCS)

# Main program sources
MAIN_SRCS := $(SRC_DIRS)/SSD_FC.cpp $(SRC_DIRS)/ldpc_codec.cpp
SRC2_MAIN_SRCS := $(SRC2_DIR)/SSD_FC.cpp $(SRC2_DIR)/ldpc_codec.cpp
TEST_MAIN_SRCS := $(SRC_DIRS)/SSD_FC.cpp $(SRC_DIRS)/ldpc_codec_test.cpp
TEST2_MAIN_SRCS := $(SRC_DIRS)/SSD_FC.cpp $(SRC_DIRS)/ldpc_codec_test2.cpp
TEST3_MAIN_SRCS := $(SRC_DIRS)/SSD_FC.cpp $(SRC_DIRS)/ldpc_codec_test3.cpp

# Object files
COMMON_C_OBJS := $(COMMON_C_SRCS:$(SRC_DIRS)/%.c=$(BUILD_DIR)/$(SRC_DIRS)/%.c.o)
COMMON_CPP_OBJS := $(COMMON_CPP_SRCS:$(SRC_DIRS)/%.cpp=$(BUILD_DIR)/$(SRC_DIRS)/%.cpp.o)
SRC2_COMMON_CPP_OBJS := $(SRC2_COMMON_CPP_SRCS:$(SRC_DIRS)/%.cpp=$(BUILD_DIR)/src2_common/%.cpp.o)
MAIN_OBJS := $(MAIN_SRCS:$(SRC_DIRS)/%.cpp=$(BUILD_DIR)/$(SRC_DIRS)/%.cpp.o)
SRC2_MAIN_OBJS := $(SRC2_MAIN_SRCS:$(SRC2_DIR)/%.cpp=$(BUILD_DIR)/$(SRC2_DIR)/%.cpp.o)
TEST_MAIN_OBJS := $(TEST_MAIN_SRCS:$(SRC_DIRS)/%.cpp=$(BUILD_DIR)/$(SRC_DIRS)/%.cpp.o)
TEST2_MAIN_OBJS := $(TEST2_MAIN_SRCS:$(SRC_DIRS)/%.cpp=$(BUILD_DIR)/$(SRC_DIRS)/%.cpp.o)
TEST3_MAIN_OBJS := $(TEST3_MAIN_SRCS:$(SRC_DIRS)/%.cpp=$(BUILD_DIR)/$(SRC_DIRS)/%.cpp.o)
OBJS := $(COMMON_C_OBJS) $(COMMON_CPP_OBJS) $(MAIN_OBJS)
SRC2_OBJS := $(COMMON_C_OBJS) $(SRC2_COMMON_CPP_OBJS) $(SRC2_MAIN_OBJS)
TEST_OBJS := $(COMMON_C_OBJS) $(COMMON_CPP_OBJS) $(TEST_MAIN_OBJS)
TEST2_OBJS := $(COMMON_C_OBJS) $(COMMON_CPP_OBJS) $(TEST2_MAIN_OBJS)
TEST3_OBJS := $(COMMON_C_OBJS) $(COMMON_CPP_OBJS) $(TEST3_MAIN_OBJS)

INC_DIRS := $(shell find $(SRC_DIRS) -type d)
INC_FLAGS := $(addprefix -I,$(INC_DIRS))

CPPFLAGS ?= $(INC_FLAGS) -MMD -MP -g -D_LDPC_DUMP
SRC2_CPPFLAGS := $(CPPFLAGS) -I$(SRC2_DIR) -DIBEX_LDPC_HEADER=\"../src2/ldpc_codec.h\"
LDFLAGS ?= -lm

# Default target: build the program
.PHONY: all
all: $(TARGET_EXEC)

# Build executable
$(TARGET_EXEC): $(OBJS)
	$(CXX) $(OBJS) -o $@ $(LDFLAGS)

ssd_fc_src2: $(SRC2_OBJS)
	$(CXX) $(SRC2_OBJS) -o $@ $(LDFLAGS)

ssd_fc_test: $(TEST_OBJS)
	$(CXX) $(TEST_OBJS) -o $@ $(LDFLAGS)

ssd_fc_test2: $(TEST2_OBJS)
	$(CXX) $(TEST2_OBJS) -o $@ $(LDFLAGS)

ssd_fc_test3: $(TEST3_OBJS)
	$(CXX) $(TEST3_OBJS) -o $@ $(LDFLAGS)

# C source compilation
$(BUILD_DIR)/$(SRC_DIRS)/%.c.o: $(SRC_DIRS)/%.c
	$(MKDIR_P) $(dir $@)
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -g -o $@

# C++ source compilation
$(BUILD_DIR)/$(SRC_DIRS)/%.cpp.o: $(SRC_DIRS)/%.cpp
	$(MKDIR_P) $(dir $@)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -g -o $@

$(BUILD_DIR)/$(SRC2_DIR)/%.cpp.o: $(SRC2_DIR)/%.cpp
	$(MKDIR_P) $(dir $@)
	$(CXX) $(SRC2_CPPFLAGS) $(CXXFLAGS) -c $< -g -o $@

$(BUILD_DIR)/src2_common/%.cpp.o: $(SRC_DIRS)/%.cpp
	$(MKDIR_P) $(dir $@)
	$(CXX) $(SRC2_CPPFLAGS) $(CXXFLAGS) -c $< -g -o $@

.PHONY: clean
clean:
	$(RM) -r $(BUILD_DIR) $(TARGET_EXEC) ssd_fc_src2
	$(RM) $(DUMPS_DIR)/*_matrix.txt
	$(RM) $(DUMPS_DIR)/*_matrix_*.txt

.PHONY: run
run:
	./$(TARGET_EXEC) LDPC config/Ibex_hd_row8.cnfg AWGN 5.30 | tee perf/test.log

MKDIR_P ?= mkdir -p

-include $(OBJS:.o=.d)
-include $(SRC2_OBJS:.o=.d)
