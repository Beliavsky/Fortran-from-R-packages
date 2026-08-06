.PHONY: test checked optimized example clean

checked:
	./tools/run_tests.sh checked

test: checked

optimized:
	./tools/run_tests.sh optimized

example: checked

clean:
	rm -rf build
