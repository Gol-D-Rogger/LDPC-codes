#Makefile
CC = g++
CXX = g++

# Two targets for two simulation programs
TARGET_EXEC ?= ssd_fc
TARGET_EXEC_DQ ?= ssd_fc_dq

BUILD_DIR ?= ./build
BUILD_DIR_DQ ?= ./build_dq
SRC_DIRS ?= ./src
DUMPS_DIR ?= ./output

# Common source files (shared by both programs)
COMMON_C_SRCS := $(shell find $(SRC_DIRS) -name '*.c')
COMMON_CPP_SRCS := $(SRC_DIRS)/fc_dsp.cpp $(SRC_DIRS)/transceiver.cpp

# Original version specific sources
ORIG_SRCS := $(SRC_DIRS)/SSD_FC.cpp $(SRC_DIRS)/ldpc_codec.cpp

# DQ version specific sources
DQ_SRCS := $(SRC_DIRS)/SSD_FC_dq.cpp

# Object files for original version
COMMON_C_OBJS := $(COMMON_C_SRCS:$(SRC_DIRS)/%.c=$(BUILD_DIR)/$(SRC_DIRS)/%.c.o)
COMMON_CPP_OBJS := $(COMMON_CPP_SRCS:$(SRC_DIRS)/%.cpp=$(BUILD_DIR)/$(SRC_DIRS)/%.cpp.o)
ORIG_OBJS := $(ORIG_SRCS:$(SRC_DIRS)/%.cpp=$(BUILD_DIR)/$(SRC_DIRS)/%.cpp.o)
OBJS := $(COMMON_C_OBJS) $(COMMON_CPP_OBJS) $(ORIG_OBJS)

# Object files for DQ version
COMMON_C_OBJS_DQ := $(COMMON_C_SRCS:$(SRC_DIRS)/%.c=$(BUILD_DIR_DQ)/$(SRC_DIRS)/%.c.o)
COMMON_CPP_OBJS_DQ := $(COMMON_CPP_SRCS:$(SRC_DIRS)/%.cpp=$(BUILD_DIR_DQ)/$(SRC_DIRS)/%.cpp.o)
DQ_OBJS := $(DQ_SRCS:$(SRC_DIRS)/%.cpp=$(BUILD_DIR_DQ)/$(SRC_DIRS)/%.cpp.o)
OBJS_DQ := $(COMMON_C_OBJS_DQ) $(COMMON_CPP_OBJS_DQ) $(DQ_OBJS) $(BUILD_DIR_DQ)/$(SRC_DIRS)/ldpc_codec_dq.cpp.o

INC_DIRS := $(shell find $(SRC_DIRS) -type d)
INC_FLAGS := $(addprefix -I,$(INC_DIRS))

CPPFLAGS ?= $(INC_FLAGS) -MMD -MP -g
CPPFLAGS_DQ ?= $(INC_FLAGS) -MMD -MP -g -DDQ_SIM

# Default target: build both programs
.PHONY: all
all: $(TARGET_EXEC) $(TARGET_EXEC_DQ)

# Build original version
$(TARGET_EXEC): $(OBJS)
	$(CXX) $(OBJS) -o $@ $(LDFLAGS) -lm

# Build DQ version
$(TARGET_EXEC_DQ): $(OBJS_DQ)
	$(CXX) $(OBJS_DQ) -o $@ $(LDFLAGS) -lm

# C source for original version
$(BUILD_DIR)/$(SRC_DIRS)/%.c.o: $(SRC_DIRS)/%.c
	$(MKDIR_P) $(dir $@)
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -g -o $@

# C++ source for original version
$(BUILD_DIR)/$(SRC_DIRS)/%.cpp.o: $(SRC_DIRS)/%.cpp
	$(MKDIR_P) $(dir $@)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -g -o $@

# C source for DQ version
$(BUILD_DIR_DQ)/$(SRC_DIRS)/%.c.o: $(SRC_DIRS)/%.c
	$(MKDIR_P) $(dir $@)
	$(CC) $(CPPFLAGS_DQ) $(CFLAGS) -c $< -g -o $@

# C++ source for DQ version
$(BUILD_DIR_DQ)/$(SRC_DIRS)/%.cpp.o: $(SRC_DIRS)/%.cpp
	$(MKDIR_P) $(dir $@)
	$(CXX) $(CPPFLAGS_DQ) $(CXXFLAGS) -c $< -g -o $@

# Special rule for ldpc_codec dq.cpp (file with space in name)
$(BUILD_DIR_DQ)/$(SRC_DIRS)/ldpc_codec_dq.cpp.o: $(SRC_DIRS)/ldpc_codec_dq.cpp
	$(MKDIR_P) $(dir $@)
	$(CXX) $(CPPFLAGS_DQ) $(CXXFLAGS) -c "$<" -g -o "$@"

.PHONY: clean
clean:
	$(RM) -r $(BUILD_DIR) $(BUILD_DIR_DQ) $(TARGET_EXEC) $(TARGET_EXEC_DQ)
	$(RM) $(DUMPS_DIR)/*_matrix.txt
	$(RM) $(DUMPS_DIR)/*_matrix_*.txt

.PHONY: clean-orig
clean-orig:
	$(RM) -r $(BUILD_DIR) $(TARGET_EXEC)

.PHONY: clean-dq
clean-dq:
	$(RM) -r $(BUILD_DIR_DQ) $(TARGET_EXEC_DQ)

.PHONY: run
run:
	./$(TARGET_EXEC) LDPC config/sdec.cnfg AWGN 5.30 | tee perf/test.log

.PHONY: run-dq
run-dq:
	./$(TARGET_EXEC_DQ) LDPC config/sdec.cnfg AWGN 5.30 | tee perf/test_dq.log

MKDIR_P ?= mkdir -p