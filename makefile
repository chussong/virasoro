UNAME := $(shell uname -s)
ifeq ($(UNAME),Linux)
	SYS = Linux
endif
ifeq ($(UNAME),Darwin)
	SYS = MacOSX-x86-64
endif
CC=g++
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
LDFLAGS=-lpthread -lgmpxx -lmpc -lmpfr -lgmp
VLDFLAGS=
ifeq ($(UNAME),Linux)
	WLDFLAGS=-L$(WSTPDIR) -lWSTP64i4 -lm -lrt -lstdc++ -ldl -luuid
endif
ifeq ($(UNAME),Darwin)
	WLDFLAGS=-L$(WSTPDIR) -lWSTPi4 -lc++ -framework Foundation
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
_DEPS=mpfc_class.h cpqmn.h hmn.h runfile.h
DEPS=$(patsubst %,$(IDIR)/%,$(_DEPS))
_VDEPS=virasoro.h
VDEPS=$(patsubst %,$(IDIR)/%,$(_VDEPS))
WDEPS=vwstp.tm
CFG=config.txt
EXECUTABLE=virasoro
WSTP=vwstp
ALLRECIPE= $(SOURCES) $(DEPS) $(EXECUTABLE) $(CFG)

ifdef MATHDIR
	ALLRECIPE += $(WSTP)
	LDFLAGS += $(WLDFLAGS)
endif

all: $(ALLRECIPE)

$(EXECUTABLE): $(OBJECTS) $(VOBJECTS) $(CFG)
	$(CC) $(OBJECTS) $(VOBJECTS) $(LDFLAGS) $(VLDFLAGS) -o $@

$(WSTP): $(OBJECTS) $(WOBJECTS)
	$(CC) $(OBJECTS) $(WOBJECTS) $(LDFLAGS) -o $@

$(ODIR)/virasoro.o: $(VSOURCES) $(DEPS) $(VDEPS) | $(ODIR)
	$(CC) $(CFLAGS) $(VCFLAGS) $(WCFLAGS) -std=c++14 $< -o $@

$(ODIR)/vwstp.o: vwstp.cpp $(DEPS) $(WDEPS) | $(ODIR)
	$(CC) $(CFLAGS) $(WCFLAGS) -std=c++14 $< -o $@

$(ODIR)/vwstptm.o: vwstptm.c | $(ODIR)
	$(CC) $(CFLAGS) $(WCFLAGS) -x c $< -o $@

vwstptm.c: $(WDEPS)
	$(WSPREP) $? -o $@

$(ODIR)/runfile.o: $(SOURCES) $(DEPS) | $(ODIR)
	$(CC) $(CFLAGS) $(WCFLAGS) -std=c++14 $< -o $@

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
	-rm -f $(EXECUTABLE) $(WSTP) $(OBJECTS) $(VOBJECTS) $(WOBJECTS) vwstptm.c
virasoro.o :
	make obj/virasoro.o
vwstp.o :
	make obj/vwstp.o
runfile.o :
	make obj/runfile.o
