
CXX = g++

CXXFLAGS = -Wall -Wno-unused -I. -Imintlib
LDFLAGS = -fPIC


UNAME = $(shell uname)

ifeq ($(UNAME), Linux)
	CXXFLAGS += -DUNIX
endif

ifeq ($(UNAME), Darwin)
	CXXFLAGS += -DUNIX -DMAC_OSX
endif

ifeq ($(debug),yes)
CXXFLAGS += -ggdb
else
CXXFLAGS += -O2
endif


#OBJECTS = mintlib/minString.o mintlib/minMap.o suffix_array.o
OBJECTS = mintlib/minString.o mintlib/minMap.o bioSequence.o

first: all

all: cd-hit-dup cd-hit-lap read-linker

.SUFFIXES: .c .obj .cpp .cc .cxx .C

.cxx.o:
	$(CXX) -c $(CXXFLAGS) -o $@ $<

cd-hit-dup: $(OBJECTS) cdhit-dup.o
	$(CXX) $(LDFLAGS) $(OBJECTS) cdhit-dup.o -o cd-hit-dup

cd-hit-lap: $(OBJECTS) cdhit-lap.o
	$(CXX) $(LDFLAGS) $(OBJECTS) cdhit-lap.o -o cd-hit-lap

read-linker: $(OBJECTS) read-linker.o
	$(CXX) $(LDFLAGS) $(OBJECTS) read-linker.o -o read-linker

clean:
	rm $(OBJECTS) cdhit-dup.o cdhit-lap.o read-linker.o
