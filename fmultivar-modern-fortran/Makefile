.PHONY: all check release-check clean

all: check

check:
	./scripts/build_and_test.sh debug

release-check:
	./scripts/build_and_test.sh release

clean:
	rm -rf build
