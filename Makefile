GENERATOR ?= Ninja
CONFIG ?= Debug

.PHONY: all build test test_failed install clean plot golden golden_run dist
all: build

build/CMakeCache.txt:
	cmake -S . -B build -G "$(GENERATOR)" -DCMAKE_BUILD_TYPE=$(CONFIG)

build: build/CMakeCache.txt
	cmake --build build

install: build
	cmake --install build

test: build
	ctest --test-dir build/test --output-on-failure -L quick

test_failed: build
	ctest --test-dir build/test --output-on-failure -V --rerun-failed

test_slow: build
	ctest --test-dir build/test --output-on-failure -L slow

test_all: build
	ctest --test-dir build/test --output-on-failure -L "quick|slow"

plot: build
	ctest --test-dir build/test --output-on-failure -L plot

external: build
	ctest --test-dir build/test -V -L external

current: build
	ctest --test-dir build/test --output-on-failure -V -L current

golden: build
	ctest --test-dir build/test --output-on-failure -V -L golden

golden_run: build
	ctest --test-dir build/test --output-on-failure -V -R GoldenRecordRun

golden_update: build/test/golden/rabe.nc
	cp $< test/golden/expected/.

clean:
	rm -rf build

dist:
	git archive --format=tar.gz --prefix=rabe/ HEAD -o rabe.tar.gz
