.PHONY: check debug-check release-check license clean

check: license debug-check

debug-check:
	./run_checks.sh debug

release-check:
	./run_checks.sh release

license:
	./test/check_license.sh

clean:
	rm -rf build
