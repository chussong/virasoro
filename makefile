CC=g++
IDIR =./include
CFLAGS=-Wall -Werror -I$(IDIR) -O3 -pg -std=c++14 -c
LDFLAGS= -lpthread -lgmpxx -lgmp
SOURCES=virasoro.cpp cpqmn.cpp hmn.cpp
ODIR = obj
_OBJECTS=$(SOURCES:.cpp=.o)
OBJECTS=$(patsubst %,$(ODIR)/%,$(_OBJECTS))
_DEPS = virasoro.h cpqmn.h hmn.h
DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))
EXECUTABLE=virasoro

all: $(SOURCES) $(DEPS) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(OBJECTS) $(LDFLAGS) -o $@

obj/virasoro.o: $(SOURCES) $(DEPS)
	$(CC) $(CFLAGS) $< -o $@

obj/cpqmn.o: cpqmn.cpp $(IDIR)/cpqmn.h
	$(CC) $(CFLAGS) $< -o $@

obj/hmn.o: hmn.cpp $(IDIR)/hmn.h cpqmn.cpp $(IDIR)/cpqmn.h
	$(CC) $(CFLAGS) $< -o $@

.PHONY: clean, virasoro.o, cpqmn.o, hmn.o
clean :
	rm $(EXECUTABLE) $(OBJECTS)
virasoro.o :
	make obj/virasoro.o
cpqmn.o :
	make obj/cpqmn.o
hmn.o :
	make obj/hmn.o
