#Makefile
CC = mpic++
CXX= mpic++

TARGET_EXEC ?= ssd_fc

BUILD_DIR ?= ./build
SRC_DIRS ?= ./src
DUMP_DIRS ?= . ./output

SRCS := $(shell find $(SRC_DIRS) -name '*.cpp' -or -name '*.c')
OBJS := $(SRCS:%=$(BUILD_DIR)/%.o)
DEPS := $(OBJS:.o=.d)

INC_DIRS := $(shell find $(SRC_DIRS) -type d)
INC_FLAGS := $(addprefix -I,$(INC_DIRS))

#CPPFLAGS ?= $(INC_FLAGS) -MMD -MP -g g -D LDPC_INIT_SYND_WT -D LDPC_DUMP
CPPFLAGS ?= $(INC_FLAGS) -MMD -MP -std=c++11
#CPPFLAGS ?= $(INC_FLAGS)-D LDPC_DUMP

$(TARGET_EXEC): $(OBJS)
		$(CXX) $(OBJS) -o $@ $(LDFLAGS)

# c source
$(BUILD_DIR)/%.c.o: %.c
	$(MKDIR_P) $(dir $@)
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -g -o $@

#c++ source
$(BUILD_DIR)/%.cpp.o: %.cpp
	$(MKDIR_P) $(dir $@)
	# $(CXX) $(CPPFLAGS) $(CXXFLAGS) - c$< -o $@
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -g -o $@
.PHONY: clean

clean:
	$(RM) -r $(BUILD_DIR)
	$(RM) $(TARGET_EXEC)
	#$(RM) $(DUMP_DIRS)/*_matrix.txt \
	#$(RM) $(DUMP_DIRS)/*_matrix_*.txt
	$(RM) $(foreach d,$(DUMP_DIRS),$(d)/*_matrix.txt $(d)/*_matrix_*.txt)

run:
	mpirun -np 1 ./ssd_fc config/Ibex_rdec_7r_4112_672.cnfg AWGN 0.02
clean_run:
	$(RM) *.log* *core*

-include $(DEPS)

MKDIR_P ?= mkdir -p
