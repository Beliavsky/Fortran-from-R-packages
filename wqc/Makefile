.PHONY: all check release clean

all: check

check:
	./build_and_test.sh

release:
	WQC_BUILD_MODE=release ./build_and_test.sh

clean:
	rm -rf build *.mod *.o
