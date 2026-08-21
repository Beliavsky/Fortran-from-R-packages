.PHONY: all debug release check clean
all: debug
debug:
	./scripts/build.sh debug
release:
	./scripts/build.sh release
check:
	./test/run_all.sh debug
clean:
	rm -rf build
