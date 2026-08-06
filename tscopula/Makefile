.PHONY: test test-check test-opt example clean

test: test-check

test-check:
	./scripts/run_tests.sh check

test-opt:
	./scripts/run_tests.sh opt

example:
	./scripts/run_example.sh opt

clean:
	rm -rf build *.o *.mod
