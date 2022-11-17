CC     ?= gcc
AR     ?= ar
PYTHON ?= python3

CFLAGS := -Wall -Wextra -g -fPIC -O2

CFLAGS += -lm

CFLAGS += -I$(shell $(PYTHON) -c "import sysconfig; print(sysconfig.get_paths()['include'])")
CFLAGS += -I$(shell $(PYTHON) -c "import numpy; print(numpy.get_include())")

ifeq ($(shell uname -s),Darwin)
CFLAGS += -undefined dynamic_lookup
endif

ifeq ($(PREFIX),)
PREFIX := /usr/local
endif

# Build C and Python libraries
.PHONY: all
all: permanent/permanent.so libpermanent.a libpermanent.so

# Install C libraries
.PHONY: install
install: permanent/tuning.h libpermanent.a libpermanent.so
	install -d $(DESTDIR)$(PREFIX)/lib/
	install -m 444 libpermanent.a $(DESTDIR)$(PREFIX)/lib/
	install -m 555 libpermanent.so $(DESTDIR)$(PREFIX)/lib/
	install -d $(DESTDIR)$(PREFIX)/include/permanent/
	install -m 444 permanent/permanent.h $(DESTDIR)$(PREFIX)/include/permanent/
	install -m 444 permanent/tuning.h $(DESTDIR)$(PREFIX)/include/permanent/

# Run tests
.PHONY: test
test: permanent/permanent.so
	$(PYTHON) -m pytest -v .

# Clean directory
.PHONY: clean
clean:
	rm -f permanent/run_tuning permanent/tuning.h
	rm -f permanent/permanent.so libpermanent.o libpermanent.a libpermanent.so
	rm -f fast_permanent.csv

# Tuning utility
permanent/run_tuning:
	$(CC) $(CFLAGS) -o $@ permanent/permanent.c permanent/run_tuning.c

# Tuning parameters
permanent/tuning.h: permanent/run_tuning
	./$^ > $@

# Python library
permanent/permanent.so: permanent/tuning.h
	$(CC) $(CFLAGS) -DTUNING_FILE=1 -shared -o $@ permanent/permanent.c permanent/py_permanent.c

# C object code
libpermanent.o: permanent/tuning.h
	$(CC) $(CFLAGS) -DTUNING_FILE=1 -c -o $@ permanent/permanent.c

# C static library
libpermanent.a: libpermanent.o
	$(AR) crs $@ $^

# C shared library
libpermanent.so: libpermanent.o
	$(CC) $(CFLAGS) -shared -o $@ $^
