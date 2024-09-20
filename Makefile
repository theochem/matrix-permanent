CXX    ?= g++
AR     ?= ar
PYTHON ?= python3

CXXFLAGS := -std=c++20 -Wall -Wextra -g -fPIC -O3

ifeq ($(shell uname -s),Darwin)
CXXFLAGS += -undefined dynamic_lookup
endif

ifneq ($(BUILD_NATIVE),)
CXXFLAGS += -march=native -mtune=native
endif

ifneq ($(RUN_TUNING),)
CXXFLAGS += -DRUN_TUNING=1
endif

# Build Python library
.PHONY: all
all: permanent/permanent.so

# Run tests
.PHONY: test
test: permanent/permanent.so
	$(PYTHON) -m pytest -v .

# Clean directory
.PHONY: clean
clean:
	rm -f src/tuning src/tuning.h src/tuning.csv permanent/permanent.so

# compile_flags.txt (clangd)
compile_flags.txt:
	echo "$(CXXFLAGS)" | sed 's/ /\n/g' > $@

# Find tuning parameters
src/tuning.h: src/permanent.h src/tuning.cc tools/tuning.py
	$(CXX) $(CXXFLAGS) -o src/tuning src/tuning.cc
	@if [ -n "$(RUN_TUNING)" ]; then \
		echo "running tuning..."; \
		src/tuning; \
		echo "writing custom tuning.h"; \
		$(PYTHON) tools/tuning.py; \
	else \
		echo "writing default tuning.h"; \
		src/tuning; \
	fi

# Compile Python library
permanent/permanent.so: src/tuning.h src/permanent.h src/py_permanent.cc tools/include_dirs.py
	$(CXX) $(CXXFLAGS) -DWITH_TUNING_FILE=1 \
		$(shell $(PYTHON) tools/include_dirs.py) \
		-shared -o $@ src/py_permanent.cc
