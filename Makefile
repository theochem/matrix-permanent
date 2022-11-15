CC     ?= gcc
PYTHON ?= python3

CFLAGS  = -fPIC -Wall -Wextra -O2 -g -lm
CFLAGS += -I$(shell $(PYTHON) -c "import sysconfig; print(sysconfig.get_paths()['include'])")
CFLAGS += -I$(shell $(PYTHON) -c "import numpy; print(numpy.get_include())")

ifeq ($(shell uname -s),Darwin)
CFLAGS += -undefined dynamic_lookup
endif

# Build C library
.PHONY: all
all: libpermanent.so

# Build Python library
.PHONY: python
python: permanent/permanent.so

# Install C library
.PHONY: install
install:
	@echo "TODO"

# Clean directory
.PHONY: clean
clean:
	rm -f libpermanent.so permanent/permanent.so permanent/run_tuning permanent/tuning.h

# Tuning utility
permanent/run_tuning:
	$(CC) $(CFLAGS) -o $@ permanent/permanent.c permanent/run_tuning.c

# Tuning parameters
permanent/tuning.h: permanent/run_tuning
	./$^ > $@

# C library
libpermanent.so: permanent/tuning.h
	$(CC) $(CFLAGS) -DTUNING_FILE=1 -shared -o $@ permanent/permanent.c

# Python library
permanent/permanent.so: permanent/tuning.h
	$(CC) $(CFLAGS) -DTUNING_FILE=1 -shared -o $@ permanent/permanent.c permanent/py_permanent.c
