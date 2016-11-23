CC=g++
IDIR =./include
CFLAGS=-Wall -Werror -Wextra -pedantic -I$(IDIR) -O3 -pg -std=c++14 -c
LDFLAGS= -lpthread -lgmpxx -lmpc -lmpfr -lgmp
SOURCES=virasoro.cpp runfile.cpp
ODIR = obj
_OBJECTS=$(SOURCES:.cpp=.o)
OBJECTS=$(patsubst %,$(ODIR)/%,$(_OBJECTS))
_DEPS = virasoro.h mpfc_class.h cpqmn.h hmn.h runfile.h
DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))
EXECUTABLE=virasoro

all: $(SOURCES) $(DEPS) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(OBJECTS) $(LDFLAGS) -o $@

obj/virasoro.o: $(SOURCES) $(DEPS) | $(ODIR)
	$(CC) $(CFLAGS) $< -o $@

obj/runfile.o: runfile.cpp $(IDIR)/runfile.h | $(ODIR)
	$(CC) $(CFLAGS) $< -o $@

$(ODIR):
	mkdir $(ODIR)

.PHONY: clean, virasoro.o
clean :
	rm $(EXECUTABLE) $(OBJECTS)
virasoro.o :
	make obj/virasoro.o
runfile.o :
	make obj/runfile.o
