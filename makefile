CC=g++
CFLAGS=-Wall -Werror -I$(IDIR) -std=c++14 -c
LDFLAGS= -lquadmath -lpthread
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
	$(CC) $(OBJECTS) $(LDFLAGS) -o $@

$(OBJECTS): $(SOURCES) $(DEPS)
	$(CC) $(CFLAGS) $< -o $@

.PHONY: clean
clean :
	rm $(OBJECTS)
