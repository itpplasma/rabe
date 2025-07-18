GENERATOR ?= Ninja
CONFIG ?= Debug

.PHONY: all build test install clean plot
all: build

build/CMakeCache.txt:
	cmake -S . -B build -G "$(GENERATOR)" -DCMAKE_BUILD_TYPE=$(CONFIG)

build: build/CMakeCache.txt
	cmake --build build

install: build
	cmake --install build

test: build
	ctest --test-dir build/test --output-on-failure -L quick

test_slow: build
	ctest --test-dir build/test --output-on-failure -L slow

test_all: build
	ctest --test-dir build/test --output-on-failure -L "quick|slow"

plot: build
	ctest --test-dir build/test --output-on-failure -L plot

external: build
	ctest --test-dir build/test -V -L external

current: build
	ctest --test-dir build/test --output-on-failure -L current

clean:
	rm -rf build
