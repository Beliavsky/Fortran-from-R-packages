.PHONY: all test optimized demo clean

all: test

test:
	./scripts/build.sh checked
	@for executable in build-checked/bin/test_*; do $$executable; done

optimized:
	./scripts/build.sh optimized
	@for executable in build-optimized/bin/test_*; do $$executable; done

demo: test
	./build-checked/bin/sharpe_rratio_demo

clean:
	rm -rf build build-checked build-optimized build-validation-*
