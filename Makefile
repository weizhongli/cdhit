
CC = g++ -Wall -ggdb
CC = g++ -pg
CC = g++

# without OpenMP

# with OpenMP
# in command line: 
# make openmp=yes
ifeq ($(openmp),no)
  CCFLAGS = -DNO_OPENMP
else
  CCFLAGS = -fopenmp
endif

# support debugging
# in command line:
# make debug=yes
# make openmp=yes debug=yes
ifeq ($(debug),yes)
CCFLAGS += -ggdb
else
CCFLAGS += -O2
endif

ifdef MAX_SEQ
CCFLAGS += -DMAX_SEQ=$(MAX_SEQ)
endif

#LDFLAGS = -static -o
LDFLAGS += -o

PROGS = cd-hit cd-hit-est cd-hit-2d cd-hit-est-2d cd-hit-div cd-hit-454

# Propagate hardening flags
CCFLAGS := $(CPPFLAGS) $(CCFLAGS) $(CXXFLAGS)

.cpp.o:
	$(CC) $(CCFLAGS) -c $<

all: create-bin-dir $(PROGS)

clean:
	rm -Rf *.o bin/

create-bin-dir: 
	mkdir -p bin

# programs

cd-hit: cdhit-common.o cdhit-utility.o cdhit.o
	$(CC) $(CCFLAGS) cdhit.o cdhit-common.o cdhit-utility.o $(LDFLAGS) bin/cd-hit

cd-hit-2d: cdhit-common.o cdhit-utility.o cdhit-2d.o
	$(CC) $(CCFLAGS) cdhit-2d.o cdhit-common.o cdhit-utility.o $(LDFLAGS) bin/cd-hit-2d

cd-hit-est: cdhit-common.o cdhit-utility.o cdhit-est.o
	$(CC) $(CCFLAGS) cdhit-est.o cdhit-common.o cdhit-utility.o $(LDFLAGS) bin/cd-hit-est

cd-hit-est-2d: cdhit-common.o cdhit-utility.o cdhit-est-2d.o
	$(CC) $(CCFLAGS) cdhit-est-2d.o cdhit-common.o cdhit-utility.o $(LDFLAGS) bin/cd-hit-est-2d

cd-hit-div: cdhit-common.o cdhit-utility.o cdhit-div.o
	$(CC) $(CCFLAGS) cdhit-div.o cdhit-common.o cdhit-utility.o $(LDFLAGS) bin/cd-hit-div

cd-hit-454: cdhit-common.o cdhit-utility.o cdhit-454.o
	$(CC) $(CCFLAGS) cdhit-454.o cdhit-common.o cdhit-utility.o $(LDFLAGS) bin/cd-hit-454

# objects
cdhit-common.o: src/cdhit-common.cpp src/cdhit-common.h
	$(CC) $(CCFLAGS) src/cdhit-common.cpp -c

cdhit-utility.o: src/cdhit-utility.cpp src/cdhit-utility.h
	$(CC) $(CCFLAGS) src/cdhit-utility.cpp -c

cdhit.o: src/cdhit.cpp src/cdhit-utility.h
	$(CC) $(CCFLAGS) src/cdhit.cpp -c

cdhit-2d.o: src/cdhit-2d.cpp src/cdhit-utility.h
	$(CC) $(CCFLAGS) src/cdhit-2d.cpp -c

cdhit-est.o: src/cdhit-est.cpp src/cdhit-utility.h
	$(CC) $(CCFLAGS) src/cdhit-est.cpp -c

cdhit-est-2d.o: src/cdhit-est-2d.cpp src/cdhit-utility.h
	$(CC) $(CCFLAGS) src/cdhit-est-2d.cpp -c

cdhit-div.o: src/cdhit-div.cpp src/cdhit-common.h
	$(CC) $(CCFLAGS) src/cdhit-div.cpp -c

cdhit-454.o: src/cdhit-454.cpp src/cdhit-common.h
	$(CC) $(CCFLAGS) src/cdhit-454.cpp -c

PREFIX ?= /usr/local/bin

install:
	for prog in $(PROGS); do \
		install -m 0755 $$prog $(PREFIX); \
	done
	install -m 0755 *.pl $(PREFIX);
