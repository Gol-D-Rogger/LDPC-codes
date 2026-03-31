CXX := g++

BUILD_DIR := build
BIN := $(BUILD_DIR)/mp_overflow_asan
CMODEL_BIN := $(BUILD_DIR)/mp_cmodel_asan
DV_LUT_BIN := $(BUILD_DIR)/dv_official_lut

SAN_FLAGS := -fsanitize=address,undefined -fno-omit-frame-pointer
COMMON_WARN := -Wall -Wextra -Wno-unused-parameter
CXXFLAGS := -std=c++17 -O1 -g $(SAN_FLAGS) $(COMMON_WARN) -I./include -I..
LDFLAGS := $(SAN_FLAGS)

MP_CODEC := ../ldpc_codec.cpp
MP_IBEX_WRAP := ../ibex_rtl_engine.cpp
MP_CMODEL := ../ldpc_c_model.c
IBEX_SRC := ../../../../IBEX/src

C_SRCS := \
	$(IBEX_SRC)/alloc.c \
	$(IBEX_SRC)/finite_lib.c \
	$(IBEX_SRC)/intio.c \
	$(IBEX_SRC)/mod2convert.c \
	$(IBEX_SRC)/mod2dense.c \
	$(IBEX_SRC)/mod2sparse.c \
	$(IBEX_SRC)/rand.c \
	$(IBEX_SRC)/vec_op.c

C_OBJS := $(patsubst %.c,$(BUILD_DIR)/%.o,$(notdir $(C_SRCS)))
COMMON_OBJS := \
	$(BUILD_DIR)/ibex_rtl_engine.o \
	$(BUILD_DIR)/ldpc_codec.o \
	$(BUILD_DIR)/transceiver.o \
	$(C_OBJS)

OVERFLOW_OBJS := $(BUILD_DIR)/mp_overflow_harness.o $(COMMON_OBJS)
CMODEL_OBJS := $(BUILD_DIR)/mp_cmodel_harness.o $(BUILD_DIR)/ldpc_c_model.o $(COMMON_OBJS)
DV_LUT_OBJS := $(BUILD_DIR)/dv_official_lut_harness.o

.PHONY: all headers asan cmodel-asan dv-official-lut clean run-config run-encode run-layer run-bf run-bf-rtl-cn run-all run-cmodel run-dv-official-lut

all: asan

headers:
	./sync_headers.sh

asan: headers $(BIN)

cmodel-asan: headers $(CMODEL_BIN)

dv-official-lut: $(DV_LUT_BIN)

$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

$(BUILD_DIR)/mp_overflow_harness.o: mp_overflow_harness.cpp | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(BUILD_DIR)/mp_cmodel_harness.o: mp_cmodel_harness.cpp | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(BUILD_DIR)/dv_official_lut_harness.o: dv_official_lut_harness.cpp | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(BUILD_DIR)/ldpc_codec.o: $(MP_CODEC) | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(BUILD_DIR)/ibex_rtl_engine.o: $(MP_IBEX_WRAP) | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(BUILD_DIR)/ldpc_c_model.o: $(MP_CMODEL) | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -include ./include/svdpi.h -x c++ -c $< -o $@

$(BUILD_DIR)/transceiver.o: $(IBEX_SRC)/transceiver.cpp | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(BUILD_DIR)/alloc.o: $(IBEX_SRC)/alloc.c | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -x c++ -c $< -o $@

$(BUILD_DIR)/finite_lib.o: $(IBEX_SRC)/finite_lib.c | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -x c++ -c $< -o $@

$(BUILD_DIR)/intio.o: $(IBEX_SRC)/intio.c | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -x c++ -c $< -o $@

$(BUILD_DIR)/mod2convert.o: $(IBEX_SRC)/mod2convert.c | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -x c++ -c $< -o $@

$(BUILD_DIR)/mod2dense.o: $(IBEX_SRC)/mod2dense.c | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -x c++ -c $< -o $@

$(BUILD_DIR)/mod2sparse.o: $(IBEX_SRC)/mod2sparse.c | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -x c++ -c $< -o $@

$(BUILD_DIR)/rand.o: $(IBEX_SRC)/rand.c | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -x c++ -c $< -o $@

$(BUILD_DIR)/vec_op.o: $(IBEX_SRC)/vec_op.c | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -x c++ -c $< -o $@

$(BIN): $(OVERFLOW_OBJS)
	$(CXX) $(OVERFLOW_OBJS) $(LDFLAGS) -o $@

$(CMODEL_BIN): $(CMODEL_OBJS)
	$(CXX) $(CMODEL_OBJS) $(LDFLAGS) -o $@

$(DV_LUT_BIN): $(DV_LUT_OBJS)
	$(CXX) $(DV_LUT_OBJS) $(LDFLAGS) -o $@

run-config: asan
	ASAN_OPTIONS=detect_leaks=0 $(BIN) --mode config

run-encode: asan
	ASAN_OPTIONS=detect_leaks=0 $(BIN) --mode encode

run-layer: asan
	ASAN_OPTIONS=detect_leaks=0 $(BIN) --mode layer

run-bf: asan
	ASAN_OPTIONS=detect_leaks=0 $(BIN) --mode bf_ibex

run-bf-rtl-cn: asan
	DVC_IBEX_FAMILY_LUT_SVH=../../../../IBEX/ibex_matrix/family_lut_output/ldpc_matrix_lut_fixed_wrap_base.svh ASAN_OPTIONS=detect_leaks=0 $(BIN) --mode bf_ibex_rtl_cn --single 10 67

run-all: asan
	ASAN_OPTIONS=detect_leaks=0 $(BIN) --mode all

run-cmodel: cmodel-asan
	ASAN_OPTIONS=detect_leaks=0 $(CMODEL_BIN)

run-dv-official-lut: dv-official-lut
	mkdir -p output
	ASAN_OPTIONS=detect_leaks=0 $(DV_LUT_BIN)

clean:
	rm -rf $(BUILD_DIR)
