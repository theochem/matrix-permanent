PREFIX ?= /usr/local

CMAKE_FLAGS :=
CMAKE_FLAGS += -DCMAKE_INSTALL_PREFIX=$(DESTDIR)$(PREFIX)
ifdef PERMANENT_TUNE
CMAKE_FLAGS += -DPERMANENT_TUNE=ON
endif
ifdef PERMANENT_PYTHON
CMAKE_FLAGS += -DPERMANENT_PYTHON=ON
CMAKE_FLAGS += -DCMAKE_EXPORT_COMPILE_COMMANDS=ON
endif

CLEAN_TARGETS :=
CLEAN_TARGETS += include/permanent/tuning.h
CLEAN_TARGETS += build dist _build _generate
CLEAN_TARGETS += permanent.*egg-info permanent.*so
CLEAN_TARGETS += compile_commands.json

.PHONY: all clean install _build

all: _build

clean:
	rm -rf $(CLEAN_TARGETS)

install: _build
	@cmake --install build

_build: build
	@cmake --build build --target all
	@if test -e build/compile_commands.json; then ln -sf build/compile_commands.json .; fi

build:
	@cmake -B build $(CMAKE_FLAGS)
