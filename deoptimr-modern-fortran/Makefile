.PHONY: debug release check clean

debug:
	./scripts/build_and_test.sh debug

release:
	./scripts/build_and_test.sh release

check: debug release

clean:
	rm -rf build
