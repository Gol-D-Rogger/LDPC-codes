#Makefile for HREModel
CC = g++
CXX = g++

# Two targets for two simulation programs
TARGET_EXEC ?= ssd_fc_hre
TARGET_EXEC_HRE_ONLY ?= ssd_fc_hre_only
TARGET_EXEC_HRE2 ?= ssd_fc_hre2

BUILD_DIR ?= ./build
BUILD_DIR_HRE_ONLY ?= ./build_hre_only
BUILD_DIR_HRE2 ?= ./build_hre2
SRC_DIRS ?= ./src
DUMPS_DIR ?= ./output

# Common source files (shared by both programs)
COMMON_C_SRCS := $(shell find $(SRC_DIRS) -name '*.c' ! -name '*.bak')
COMMON_CPP_SRCS := $(SRC_DIRS)/fc_dsp.cpp $(SRC_DIRS)/transceiver.cpp $(SRC_DIRS)/isoc23_compat.cpp

# Original version specific sources (full simulation)
ORIG_SRCS := $(SRC_DIRS)/SSD_FC.cpp $(SRC_DIRS)/ldpc_codec.cpp

# HRE-only version specific sources
HRE_ONLY_SRCS := $(SRC_DIRS)/SSD_FC_hre_only.cpp $(SRC_DIRS)/transceiver_hre_only.cpp
HRE2_MAIN_SRCS := $(SRC_DIRS)/SSD_FC_hre2.cpp $(SRC_DIRS)/ldpc_codec.cpp

# Object files for original version
COMMON_C_OBJS := $(COMMON_C_SRCS:$(SRC_DIRS)/%.c=$(BUILD_DIR)/$(SRC_DIRS)/%.c.o)
COMMON_CPP_OBJS := $(COMMON_CPP_SRCS:$(SRC_DIRS)/%.cpp=$(BUILD_DIR)/$(SRC_DIRS)/%.cpp.o)
ORIG_OBJS := $(ORIG_SRCS:$(SRC_DIRS)/%.cpp=$(BUILD_DIR)/$(SRC_DIRS)/%.cpp.o)
OBJS := $(COMMON_C_OBJS) $(COMMON_CPP_OBJS) $(ORIG_OBJS)

# Object files for HRE-only version (only need transceiver_hre_only, not full ldpc)
HRE_ONLY_C_OBJS := $(COMMON_C_SRCS:$(SRC_DIRS)/%.c=$(BUILD_DIR_HRE_ONLY)/$(SRC_DIRS)/%.c.o)
HRE_ONLY_CPP_OBJS := $(SRC_DIRS)/isoc23_compat.cpp
HRE_ONLY_CPP_OBJS := $(HRE_ONLY_CPP_OBJS:$(SRC_DIRS)/%.cpp=$(BUILD_DIR_HRE_ONLY)/$(SRC_DIRS)/%.cpp.o)
HRE_ONLY_MAIN_OBJS := $(HRE_ONLY_SRCS:$(SRC_DIRS)/%.cpp=$(BUILD_DIR_HRE_ONLY)/$(SRC_DIRS)/%.cpp.o)
OBJS_HRE_ONLY := $(HRE_ONLY_C_OBJS) $(HRE_ONLY_CPP_OBJS) $(HRE_ONLY_MAIN_OBJS)

# Object files for HRE2 version (channel only, new HRE logic)
HRE2_C_OBJS := $(COMMON_C_SRCS:$(SRC_DIRS)/%.c=$(BUILD_DIR_HRE2)/$(SRC_DIRS)/%.c.o)
HRE2_CPP_OBJS := $(COMMON_CPP_SRCS:$(SRC_DIRS)/%.cpp=$(BUILD_DIR_HRE2)/$(SRC_DIRS)/%.cpp.o)
HRE2_MAIN_OBJS := $(HRE2_MAIN_SRCS:$(SRC_DIRS)/%.cpp=$(BUILD_DIR_HRE2)/$(SRC_DIRS)/%.cpp.o)
OBJS_HRE2 := $(HRE2_C_OBJS) $(HRE2_CPP_OBJS) $(HRE2_MAIN_OBJS)

INC_DIRS := $(shell find $(SRC_DIRS) -type d)
INC_FLAGS := $(addprefix -I,$(INC_DIRS))

CPPFLAGS ?= $(INC_FLAGS) -MMD -MP -g
CPPFLAGS_HRE2 := $(CPPFLAGS) -DHRE2_MODE
LDFLAGS ?= -lm

# Default target: build both programs
.PHONY: all
all: $(TARGET_EXEC) $(TARGET_EXEC_HRE_ONLY) $(TARGET_EXEC_HRE2)
	@:

# Build original version (full HRE simulation with LDPC)
$(TARGET_EXEC): $(OBJS)
	$(CXX) $(OBJS) -o $@ $(LDFLAGS)

# Build HRE-only version (channel simulation only, no LDPC)
$(TARGET_EXEC_HRE_ONLY): $(OBJS_HRE_ONLY)
	$(CXX) $(OBJS_HRE_ONLY) -o $@ $(LDFLAGS)

# Build HRE2 version (channel simulation only, new HRE flip + optional AWGN)
$(TARGET_EXEC_HRE2): $(OBJS_HRE2)
	$(CXX) $(OBJS_HRE2) -o $@ $(LDFLAGS)

# C source for original version
$(BUILD_DIR)/$(SRC_DIRS)/%.c.o: $(SRC_DIRS)/%.c
	$(MKDIR_P) $(dir $@)
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -g -o $@

# C++ source for original version
$(BUILD_DIR)/$(SRC_DIRS)/%.cpp.o: $(SRC_DIRS)/%.cpp
	$(MKDIR_P) $(dir $@)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -g -o $@

# C source for HRE-only version
$(BUILD_DIR_HRE_ONLY)/$(SRC_DIRS)/%.c.o: $(SRC_DIRS)/%.c
	$(MKDIR_P) $(dir $@)
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -g -o $@

# C++ source for HRE-only version
$(BUILD_DIR_HRE_ONLY)/$(SRC_DIRS)/%.cpp.o: $(SRC_DIRS)/%.cpp
	$(MKDIR_P) $(dir $@)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -g -o $@

# C source for HRE2 version
$(BUILD_DIR_HRE2)/$(SRC_DIRS)/%.c.o: $(SRC_DIRS)/%.c
	$(MKDIR_P) $(dir $@)
	$(CC) $(CPPFLAGS_HRE2) $(CFLAGS) -c $< -g -o $@

# C++ source for HRE2 version
$(BUILD_DIR_HRE2)/$(SRC_DIRS)/%.cpp.o: $(SRC_DIRS)/%.cpp
	$(MKDIR_P) $(dir $@)
	$(CXX) $(CPPFLAGS_HRE2) $(CXXFLAGS) -c $< -g -o $@

.PHONY: clean
clean:
	$(RM) -r $(BUILD_DIR) $(BUILD_DIR_HRE_ONLY) $(BUILD_DIR_HRE2) $(TARGET_EXEC) $(TARGET_EXEC_HRE_ONLY) $(TARGET_EXEC_HRE2)
	$(RM) $(DUMPS_DIR)/*_matrix.txt
	$(RM) $(DUMPS_DIR)/*_matrix_*.txt

.PHONY: clean-orig
clean-orig:
	$(RM) -r $(BUILD_DIR) $(TARGET_EXEC)

.PHONY: clean-hre-only
clean-hre-only:
	$(RM) -r $(BUILD_DIR_HRE_ONLY) $(TARGET_EXEC_HRE_ONLY)

.PHONY: clean-hre2
clean-hre2:
	$(RM) -r $(BUILD_DIR_HRE2) $(TARGET_EXEC_HRE2)

.PHONY: run
run:
	./$(TARGET_EXEC) LDPC config/21x150_hre0_dis.cnfg AWGN 5.30 0 ./matrix | tee perf/test.log

.PHONY: run-hre-only
run-hre-only:
	./$(TARGET_EXEC_HRE_ONLY) LDPC config/21x150_hre0_dis.cnfg AWGN 5.30 | tee perf/test_hre_only.log

.PHONY: run-hre2
run-hre2:
	./$(TARGET_EXEC_HRE2) LDPC config/21x150_hre0_dis.cnfg AWGN 5.30 | tee perf/test_hre2.log

MKDIR_P ?= mkdir -p

-include $(OBJS:.o=.d)
-include $(OBJS_HRE_ONLY:.o=.d)
-include $(OBJS_HRE2:.o=.d)
