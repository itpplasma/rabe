CONFIG ?= Debug

# Prevent ambient shell values from silently changing the fetched libneo ref.
# Pass them explicitly: make test LIBNEO_REF=main  or  LIBNEO_PATH=/path/to/libneo
unexport LIBNEO_REF LIBNEO_PATH

# Forward only values that were explicitly supplied on the make command line.
_LIBNEO_CMAKE_ARGS :=
ifdef LIBNEO_REF
  _LIBNEO_CMAKE_ARGS += -DLIBNEO_REF=$(LIBNEO_REF)
endif
ifdef LIBNEO_PATH
  _LIBNEO_CMAKE_ARGS += -DLIBNEO_PATH=$(LIBNEO_PATH)
endif

.PHONY: all build test test_failed install clean plot golden golden_run dist internal python
all: build

build/CMakeCache.txt:
	cmake -S . -B build -DCMAKE_BUILD_TYPE=$(CONFIG) $(_LIBNEO_CMAKE_ARGS)

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

python:
	cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DBUILD_PYTHON_BINDINGS=ON
	cmake --build build

# Private test data and external comparisons (not available in public repo)
internal: build
	ctest --test-dir build/test -V -L internal
