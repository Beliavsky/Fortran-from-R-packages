.PHONY: debug release test check clean

debug:
	./scripts/build.sh debug

release:
	./scripts/build.sh release

test: debug
	./scripts/run_tests.sh debug

check:
	./test/check_license.sh
	./scripts/build.sh debug
	./scripts/run_tests.sh debug
	./scripts/build.sh release
	./scripts/run_tests.sh release

clean:
	rm -rf build
