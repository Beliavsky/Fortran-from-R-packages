FC = gfortran
AR = ar

SRC = \
  src/waveslim_kinds.f90 \
  src/waveslim_status.f90 \
  src/waveslim_types.f90 \
  src/waveslim_math.f90 \
  src/waveslim_linalg.f90 \
  src/waveslim_filters.f90 \
  src/waveslim_transform_1d.f90 \
  src/waveslim_packet.f90 \
  src/waveslim_transform_nd.f90 \
  src/waveslim_dualtree.f90 \
  src/waveslim_hilbert_stats.f90 \
  src/waveslim_statistics.f90 \
  src/waveslim_denoise.f90 \
  src/waveslim_long_memory.f90 \
  src/waveslim_extended.f90 \
  src/waveslim.f90

TEST_NAMES = \
  test_transforms_1d \
  test_packets \
  test_multidimensional \
  test_statistics_denoise \
  test_long_memory \
  test_extended \
  test_errors

CHECK_FLAGS = -std=f2018 -O0 -g -Wall -Wextra -Werror -fcheck=all \
  -fbacktrace -ffpe-trap=invalid,zero,overflow
RELEASE_FLAGS = -std=f2018 -O3 -Wall -Wextra -Werror

.PHONY: all check release app example clean
all: check

check: build/check/.tested
release: build/release/.tested
app: build/check/app/waveslim_demo
example: build/check/example/image_packet_example

build/check/.library: $(SRC)
	@mkdir -p build/check/mod build/check/obj build/check/test build/check/app build/check/example
	@cd build/check/obj && $(FC) $(CHECK_FLAGS) -J ../mod -I ../mod -c $(addprefix ../../../,$(SRC))
	@$(AR) rcs build/check/libwaveslim.a build/check/obj/*.o
	@touch $@

build/release/.library: $(SRC)
	@mkdir -p build/release/mod build/release/obj build/release/test build/release/app build/release/example
	@cd build/release/obj && $(FC) $(RELEASE_FLAGS) -J ../mod -I ../mod -c $(addprefix ../../../,$(SRC))
	@$(AR) rcs build/release/libwaveslim.a build/release/obj/*.o
	@touch $@

build/check/test_support.o: test/test_support.f90 build/check/.library
	@$(FC) $(CHECK_FLAGS) -J build/check/mod -I build/check/mod -c $< -o $@

build/release/test_support.o: test/test_support.f90 build/release/.library
	@$(FC) $(RELEASE_FLAGS) -J build/release/mod -I build/release/mod -c $< -o $@

define CHECK_TEST_template
build/check/test/$(1): test/$(1).f90 build/check/.library build/check/test_support.o
	@$$(FC) $$(CHECK_FLAGS) -J build/check/mod -I build/check/mod $$< \
	  build/check/test_support.o build/check/libwaveslim.a -o $$@
endef
$(foreach test,$(TEST_NAMES),$(eval $(call CHECK_TEST_template,$(test))))

define RELEASE_TEST_template
build/release/test/$(1): test/$(1).f90 build/release/.library build/release/test_support.o
	@$$(FC) $$(RELEASE_FLAGS) -J build/release/mod -I build/release/mod $$< \
	  build/release/test_support.o build/release/libwaveslim.a -o $$@
endef
$(foreach test,$(TEST_NAMES),$(eval $(call RELEASE_TEST_template,$(test))))

build/check/.tested: $(addprefix build/check/test/,$(TEST_NAMES))
	@for test in $^; do ./$$test || exit 1; done
	@touch $@

build/release/.tested: $(addprefix build/release/test/,$(TEST_NAMES))
	@for test in $^; do ./$$test || exit 1; done
	@touch $@

build/check/app/waveslim_demo: app/waveslim_demo.f90 build/check/.library
	@$(FC) $(CHECK_FLAGS) -J build/check/mod -I build/check/mod $< \
	  build/check/libwaveslim.a -o $@

build/check/example/image_packet_example: example/image_packet_example.f90 build/check/.library
	@$(FC) $(CHECK_FLAGS) -J build/check/mod -I build/check/mod $< \
	  build/check/libwaveslim.a -o $@

clean:
	rm -rf build *.o *.mod
