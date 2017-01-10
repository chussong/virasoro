UNAME := $(shell uname -s)
ifeq ($(UNAME),Linux)
	SYS = Linux
	CC = gcc
	CXX = g++ -std=c++14
endif
ifeq ($(UNAME),Darwin)
	SYS = MacOSX-x86-64
	CC = clang
	CXX = clang++ -std=c++14
endif
ifeq ($(UNAME),Linux)
	MATHDIR = $(shell echo $$MATHEMATICA_HOME)
endif
ifeq ($(UNAME),Darwin)
ifneq ($(wildcard /Applications/Mathematica.app/Contents/.*),)
	MATHDIR = /Applications/Mathematica.app/Contents
endif
endif
WSTPDIR = $(MATHDIR)/SystemFiles/Links/WSTP/DeveloperKit/$(SYS)/CompilerAdditions
WSPREP = $(WSTPDIR)/wsprep
IDIR =./include
CFLAGS=-Wall -Wextra -pedantic -Werror -I$(IDIR) -O3 -pg -c
VCFLAGS=
WCFLAGS=-Wno-unused-parameter -I$(WSTPDIR)
ifdef MATHDIR
	CFLAGS += -DWSTP
endif
GMP=/usr/local/lib/libgmp.a
MPFR=/usr/local/lib/libmpfr.a
MPC=/usr/local/lib/libmpc.a
LDFLAGS=-lpthread $(MPC) $(MPFR) $(GMP)
VLDFLAGS=
ifeq ($(UNAME),Darwin)
	CFLAGS+=-I/usr/local/include
	LDFLAGS+=-L/usr/local/lib
endif
ifdef MATHDIR
ifeq ($(UNAME),Linux)
	WSTP=$(WSTPDIR)/libWSTP64i4.a
	WLDFLAGS=$(WSTP) -lm -lrt -lstdc++ -ldl -luuid -no-pie
endif
ifeq ($(UNAME),Darwin)
	WSTP=$(WSTPDIR)/libWSTPi4.a
	WLDFLAGS=$(WSTP) -lc++ -framework Foundation
endif
endif
SOURCES=runfile.cpp
VSOURCES=virasoro.cpp
WSOURCES=vwstp.cpp vwstptm.c
ODIR=obj
_OBJECTS=$(SOURCES:.cpp=.o)
OBJECTS=$(patsubst %,$(ODIR)/%,$(_OBJECTS))
_VOBJECTS=$(VSOURCES:.cpp=.o)
VOBJECTS=$(patsubst %,$(ODIR)/%,$(_VOBJECTS))
_WOBJECTS=vwstp.o vwstptm.o
WOBJECTS=$(patsubst %,$(ODIR)/%,$(_WOBJECTS))
_DEPS=mpreal.h mpcomplex.h cpqmn.h hmn.h runfile.h
DEPS=$(patsubst %,$(IDIR)/%,$(_DEPS))
_VDEPS=virasoro.h
VDEPS=$(patsubst %,$(IDIR)/%,$(_VDEPS))
WDEPS=vwstp.tm
CFG=config.txt
EXECUTABLE=virasoro
VWSTP=vwstp
ALLRECIPE= $(SOURCES) $(DEPS) $(EXECUTABLE) $(CFG)

ifdef MATHDIR
	ALLRECIPE += $(VWSTP)
	LDFLAGS += $(WLDFLAGS)
endif

all: $(ALLRECIPE)

$(EXECUTABLE): $(OBJECTS) $(VOBJECTS) $(CFG)
	$(CXX) $(OBJECTS) $(VOBJECTS) $(LDFLAGS) $(VLDFLAGS) -o $@

$(VWSTP): $(OBJECTS) $(WOBJECTS)
	$(CXX) $(OBJECTS) $(WOBJECTS) $(LDFLAGS) -o $@

$(ODIR)/virasoro.o: $(VSOURCES) $(DEPS) $(VDEPS) | $(ODIR)
	$(CXX) $(CFLAGS) $(VCFLAGS) $(WCFLAGS) $< -o $@

$(ODIR)/vwstp.o: vwstp.cpp $(DEPS) $(WDEPS) | $(ODIR)
	$(CXX) $(CFLAGS) $(WCFLAGS) $< -o $@

$(ODIR)/vwstptm.o: vwstptm.c | $(ODIR)
	$(CC) $(CFLAGS) $(WCFLAGS) -x c $< -o $@

vwstptm.c: $(WDEPS)
	$(WSPREP) $? -o $@

$(ODIR)/runfile.o: $(SOURCES) $(DEPS) | $(ODIR)
	$(CXX) $(CFLAGS) $(WCFLAGS) $< -o $@

$(ODIR):
	mkdir $(ODIR)

$(CFG): virasoro.cpp
	$(file > $(CFG),[default parameters])
	$(file >> $(CFG),maxThreads=8)
	$(file >> $(CFG),precision=768)
	$(file >> $(CFG),tolerance=1e-10)
	$(file >> $(CFG),showProgressBar=true)

.PHONY: clean, virasoro.o, vwstp.o, runfile.o
clean :
	-rm -f $(EXECUTABLE) $(VWSTP) $(OBJECTS) $(VOBJECTS) $(WOBJECTS) vwstptm.c
virasoro.o :
	make obj/virasoro.o
vwstp.o :
	make obj/vwstp.o
runfile.o :
	make obj/runfile.o
