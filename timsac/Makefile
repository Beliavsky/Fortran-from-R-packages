FC = gfortran
AR = ar
FFLAGS ?= -std=f2018 -O2 -ffree-line-length-none -w
DBGFLAGS ?= -std=f2018 -O0 -g -ffree-line-length-none -Wall -Wextra -fcheck=all

BUILD := build
MODDIR := $(BUILD)/mod
LIB := $(BUILD)/libtimsac.a

LEGACY_SRC := $(wildcard src/legacy/*.f90)
LEGACY_OBJ := $(patsubst src/legacy/%.f90,$(BUILD)/%.o,$(LEGACY_SRC))
MODULE_OBJ := $(BUILD)/timsac_kinds.o $(BUILD)/timsac_raw.o $(BUILD)/timsac.o
OBJECTS := $(MODULE_OBJ) $(LEGACY_OBJ)

.PHONY: all test example debug clean

all: $(LIB)

$(BUILD) $(MODDIR):
	mkdir -p $@

$(BUILD)/timsac_kinds.o: src/timsac_kinds.f90 | $(BUILD) $(MODDIR)
	$(FC) $(FFLAGS) -J$(MODDIR) -c $< -o $@

$(BUILD)/timsac_raw.o: src/timsac_raw.f90 $(BUILD)/timsac_kinds.o | $(BUILD) $(MODDIR)
	$(FC) $(FFLAGS) -I$(MODDIR) -J$(MODDIR) -c $< -o $@

$(BUILD)/timsac.o: src/timsac.f90 $(BUILD)/timsac_kinds.o $(BUILD)/timsac_raw.o | $(BUILD) $(MODDIR)
	$(FC) $(FFLAGS) -I$(MODDIR) -J$(MODDIR) -c $< -o $@

$(BUILD)/%.o: src/legacy/%.f90 $(BUILD)/timsac_kinds.o | $(BUILD) $(MODDIR)
	$(FC) $(FFLAGS) -I$(MODDIR) -c $< -o $@

$(LIB): $(OBJECTS)
	$(AR) rcs $@ $^

test: $(LIB)
	$(FC) $(FFLAGS) -I$(MODDIR) test/test_timsac.f90 $(LIB) -o $(BUILD)/test_timsac
	$(BUILD)/test_timsac

example: $(LIB)
	$(FC) $(FFLAGS) -I$(MODDIR) app/example_autocorrelation.f90 $(LIB) -o $(BUILD)/example_autocorrelation
	$(BUILD)/example_autocorrelation

debug: FFLAGS := $(DBGFLAGS)
debug: clean test

clean:
	rm -rf $(BUILD)
