CC=g++
CFLAGS=-Wall -Werror -pthread -libquadmath -I$(IDIR) -std=c++14 -c
LDFLAGS=
SOURCES=virasoro.cpp
ODIR = obj
_OBJECTS=$(SOURCES:.cpp=.o)
OBJECTS=$(patsubst %,$(ODIR)/%,$(_OBJECTS))
IDIR =./include
_DEPS = virasoro.h
DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))
EXECUTABLE=virasoro

all: $(SOURCES) $(DEPS) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@

$(OBJECTS): $(SOURCES) $(DEPS)
	$(CC) $(CFLAGS) $< -o $@

.PHONY: clean
clean :
	rm $(OBJECTS)
