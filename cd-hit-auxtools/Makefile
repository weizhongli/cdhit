
CC = g++

CFLAGS = -Wall -Wno-unused -I. -Imintlib
LFLAGS = -fPIC


UNAME = $(shell uname)

ifeq ($(UNAME), Linux)
	CFLAGS += -DUNIX
endif

ifeq ($(UNAME), Darwin)
	CFLAGS += -DUNIX -DMAC_OSX
endif

ifeq ($(debug),yes)
CFLAGS += -ggdb
else
CFLAGS += -O2
endif


#OBJECTS = mintlib/minString.o mintlib/minMap.o suffix_array.o
OBJECTS = mintlib/minString.o mintlib/minMap.o bioSequence.o

first: all

all: cd-hit-dup cd-hit-lap read-linker

.SUFFIXES: .c .obj .cpp .cc .cxx .C

.cxx.o:
	$(CC) -c $(CFLAGS) -o $@ $<

cd-hit-dup: $(OBJECTS) cdhit-dup.o
	$(CC) $(LFLAGS) $(OBJECTS) cdhit-dup.o -o cd-hit-dup

cd-hit-lap: $(OBJECTS) cdhit-lap.o
	$(CC) $(LFLAGS) $(OBJECTS) cdhit-lap.o -o cd-hit-lap

read-linker: $(OBJECTS) read-linker.o
	$(CC) $(LFLAGS) $(OBJECTS) read-linker.o -o read-linker

clean:
	rm $(OBJECTS) cdhit-dup.o cdhit-lap.o read-linker.o
