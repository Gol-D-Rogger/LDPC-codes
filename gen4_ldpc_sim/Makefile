#Makefile
CC = g++
CXX = g++

TARGET_EXEC ?= ssd_fc

BUILD_DIR ?= ./build
SRC_DIRS ?= ./src
DUMPS_DIR ?= ./output

SRCS := $(shell find $(SRC_DIRS) -name '*.c' -or -name '*.cpp')
OBJS := $(SRCS:%=$(BUILD_DIR)/%.o)
DEPS := $(OBJS:.o=.d)

INC_DIRS := $(shell find $(SRC_DIRS) -type d)
INC_FLAGS := $(addprefix -I,$(INC_DIRS))

CPPFLAGS ?= $(INC_FLAGS) -MMD -MP -g

$(TARGET_EXEC): $(OBJS)
	$(CXX) $(OBJS) -o $@ $(LDFLAGS) -lm

# C source
$(BUILD_DIR)/%.c.o: %.c
	$(MKDIR_P) $(dir $@)
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -g -o $@

# C++ source
$(BUILD_DIR)/%.cpp.o: %.cpp
	$(MKDIR_P) $(dir $@)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -g -o $@

.PHONY: clean
clean:
	$(RM) -r $(BUILD_DIR) $(TARGET_EXEC)
	$(RM) $(DUMPS_DIR)/*_matrix.txt
	$(RM) $(DUMPS_DIR)/*_matrix_*.txt

run:
	./ssd_fc LDPC config/sdec.cnfg AWGN 5.30 | tee perf/test.log

-include $(DEPS)

MKDIR_P ?= mkdir -p