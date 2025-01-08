GENERATOR ?= Ninja

.PHONY: all build test install clean
all: build

build/CMakeCache.txt:
	cmake -S . -B build -G $(GENERATOR)

build: build/CMakeCache.txt
	cmake --build build

install: build
	cmake --install build

test: build
	ctest --test-dir build/tests --output-on-failure

clean:
	rm -rf build
