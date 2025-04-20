GENERATOR ?= Ninja
CONFIG ?= Debug

.PHONY: all build test install clean
all: build

build/CMakeCache.txt:
	cmake -S . -B build -G "$(GENERATOR)" -DCMAKE_BUILD_TYPE=$(CONFIG)

build: build/CMakeCache.txt
	cmake --build build

install: build
	cmake --install build

test: build
	ctest --test-dir build/test --output-on-failure

clean:
	rm -rf build
