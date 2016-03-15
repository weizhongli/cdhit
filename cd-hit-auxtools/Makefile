CXX ?= g++

CXXFLAGS ?= -Wall -Wno-unused
LFLAGS = -fPIC
CPPFLAGS = -I. -Imintlib

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
CXXFLAGS ?= -O2
endif


#OBJECTS = mintlib/minString.o mintlib/minMap.o suffix_array.o
OBJECTS = mintlib/minString.o mintlib/minMap.o bioSequence.o

first: all

all: cd-hit-dup cd-hit-lap read-linker

.SUFFIXES: .c .obj .cpp .cc .cxx .C

.cxx.o:
	$(CXX) -c $(CXXFLAGS) $(CPPFLAGS) -o $@ $<

cd-hit-dup: $(OBJECTS) cdhit-dup.o
	$(CXX) $(CXXFLAGS) $(LFLAGS) $(OBJECTS) cdhit-dup.o -o cd-hit-dup

cd-hit-lap: $(OBJECTS) cdhit-lap.o
	$(CXX) $(CXXFLAGS) $(LFLAGS) $(OBJECTS) cdhit-lap.o -o cd-hit-lap

read-linker: $(OBJECTS) read-linker.o
	$(CXX) $(CXXFLAGS) $(LFLAGS) $(OBJECTS) read-linker.o -o read-linker

clean:
	rm $(OBJECTS) cdhit-dup.o cdhit-lap.o read-linker.o
