.PHONY: debug release check clean

debug:
	./scripts/build_and_test.sh debug

release:
	./scripts/build_and_test.sh release

check: debug

clean:
	rm -rf build *.o *.mod
