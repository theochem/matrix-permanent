PREFIX ?= /usr/local

CMAKE_FLAGS :=
CMAKE_FLAGS += -DCMAKE_INSTALL_PREFIX=$(DESTDIR)$(PREFIX)
ifdef PERMANENT_TUNE
CMAKE_FLAGS += -DPERMANENT_TUNE=ON
endif

CLEAN_TARGETS :=
CLEAN_TARGETS += include/permanent/tuning.h
CLEAN_TARGETS += build dist _build _generate
CLEAN_TARGETS += permanent.*egg-info permanent.*so

.PHONY: all clean install _build

all: _build

clean:
	rm -rf $(CLEAN_TARGETS)

install: _build
	cmake --install build

_build: build
	cmake --build build --target all

build:
	cmake -B build $(CMAKE_FLAGS) .
