CXX    ?= g++
AR     ?= ar
PYTHON ?= python3

CXXFLAGS := -std=c++17 -Wall -Wextra -g -fPIC -O3

ifeq ($(shell uname -s),Darwin)
CXXFLAGS += -undefined dynamic_lookup
endif

ifneq ($(RUN_TUNING),)
CXXFLAGS += -DRUN_TUNING=1
endif

ifeq ($(PREFIX),)
PREFIX := /usr/local
endif

# Build C and Python libraries
.PHONY: all
all: libpermanent.a libpermanent.so permanent/permanent.so 

# Build C libraries
.PHONY: c
c: libpermanent.a libpermanent.so

# Build Python libraries
.PHONY: python
python: permanent/permanent.so

# Install C libraries
.PHONY: install
install: src/tuning.h libpermanent.a libpermanent.so
	install -d $(DESTDIR)$(PREFIX)/lib/
	install -m 444 libpermanent.a $(DESTDIR)$(PREFIX)/lib/
	install -m 555 libpermanent.so $(DESTDIR)$(PREFIX)/lib/
	install -d $(DESTDIR)$(PREFIX)/include/src/
	install -m 444 src/permanent.h $(DESTDIR)$(PREFIX)/include/src/
	install -m 444 src/tuning.h $(DESTDIR)$(PREFIX)/include/src/

# Run tests
.PHONY: test
test: permanent/permanent.so
	$(PYTHON) -m pytest -v .

# Clean directory
.PHONY: clean
clean:
	rm -f src/tuning src/tuning.h src/tuning.csv
	rm -f src/libpermanent.o permanent/permanent.so libpermanent.a libpermanent.so

# compile_flags.txt (clangd)
compile_flags.txt:
	echo "$(CXXFLAGS)" | sed 's/ /\n/g' > $@

# Find tuning parameters
src/tuning.h: src/permanent.h src/tuning.cc src/tuning.py
	$(CXX) $(CXXFLAGS) -o src/tuning src/permanent.cc src/tuning.cc
	src/tuning
	[ -n "$(RUN_TUNING)" ] && $(PYTHON) src/tuning.py || true

# Compile Python library
permanent/permanent.so: src/tuning.h src/permanent.h src/permanent.cc src/py_permanent.cc
	$(CXX) $(CXXFLAGS) -DWITH_TUNING_FILE=1 \
		-I$(shell $(PYTHON) -c "import sysconfig; print(sysconfig.get_paths()['include'])") \
		-I$(shell $(PYTHON) -c "import numpy; print(numpy.get_include())") \
		-shared -o $@ src/permanent.cc src/py_permanent.cc

# Compile object code
src/libpermanent.o: src/tuning.h src/permanent.h src/permanent.cc
	$(CXX) $(CXXFLAGS) -DWITH_TUNING_FILE=1 -c -o $@ src/permanent.cc

# Compile static library
libpermanent.a: src/libpermanent.o
	$(AR) crs $@ $^

# Compile shared library
libpermanent.so: src/libpermanent.o
	$(CXX) $(CXXFLAGS) -shared -o $@ $^


