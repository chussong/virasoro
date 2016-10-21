CC=g++
CFLAGS=-Wall -Werror -I$(IDIR) -fext-numeric-literals -O3 -std=c++14 -c
LDFLAGS= -lquadmath -lpthread -lgmp
SOURCES=virasoro128.cpp
ODIR = obj
_OBJECTS=$(SOURCES:.cpp=.o)
OBJECTS=$(patsubst %,$(ODIR)/%,$(_OBJECTS))
IDIR =./include
_DEPS = virasoro128.h
DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))
EXECUTABLE=virasoro128

all: $(SOURCES) $(DEPS) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(OBJECTS) $(LDFLAGS) -o $@

$(OBJECTS): $(SOURCES) $(DEPS)
	$(CC) $(CFLAGS) $< -o $@

.PHONY: clean
clean :
	rm $(OBJECTS)
