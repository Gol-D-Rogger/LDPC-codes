#Makefile
CC = gcc
CXX = gcc

TARGET_EXEC ?= ldpc_gen

BUILD_DIR ?= ./build
SRC_DIRS ?= ./src
DUMPS_DIR ?= ./output

SRCS := $(shell find $(SRC_DIRS) -name '*.cpp' -or -name '*.c')
OBJS := $(SRCS:%=$(BUILD_DIR)/%.o)
DEPS := $(OBJS:.o=.d)

INC_DIRS := $(shell find $(SRC_DIRS) -type d)
INC_FLAGS := $(addprefix -I,$(INC_DIRS))
#CFLAGS += -Wall -Wextra -Wmisleading-indentation -Wempty-body -Wparentheses

CPPFLAGS ?= $(INC_FLAGS) -MMD -MP -g

$(TARGET_EXEC): $(OBJS)
	$(CXX) $(OBJS) -o $@ $(LDFLAGS) -lm

# C source
$(BUILD_DIR)/%.c.o: %.c
	$(MKDIR_P) $(dir $@)
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -g -o $@

.PHONY: clean
clean:
	$(RM) -r $(BUILD_DIR) \
	$(RM) $(TARGET_EXEC)

run:
	./ldpc_gen 149 129 1
	./ldpc_gen 149 129 2
	./ldpc_gen 149 129 3
	./ldpc_gen 149 129 4
	./ldpc_gen 149 129 5

-include $(DEPS)

MKDIR_P ?= mkdir -p