
CC = g++ -Wall -ggdb
CC = g++ -pg
CC = g++

# without OpenMP
CCFLAGS = -DNO_OPENMP

# with OpenMP
# in command line: 
# make openmp=yes
ifeq ($(openmp),yes)
CCFLAGS = -fopenmp
endif

# support debugging
# in command line:
# make debug=yes
# make openmp=yes debug=yes
ifeq ($(debug),yes)
CCFLAGS += -ggdb -O0
else
CCFLAGS += -O2
endif

#LDFLAGS = -static -o
LDFLAGS = -o

PROGS = cd-hit cd-hit-est cd-hit-2d cd-hit-est-2d cd-hit-div cd-hit-454

COMMON_OBJS = cdhit-common.o cdhit-utility.o MurmurHash3.o

.c++.o:
	$(CC) $(CCFLAGS) -c $<

all: $(PROGS)

clean:
	rm *.o $(PROGS)

# programs

cd-hit: $(COMMON_OBJS) cdhit.o
	$(CC) $(CCFLAGS) cdhit.o $(COMMON_OBJS) $(LDFLAGS) cd-hit

cd-hit-2d: $(COMMON_OBJS) cdhit-2d.o
	$(CC) $(CCFLAGS) cdhit-2d.o $(COMMON_OBJS) $(LDFLAGS) cd-hit-2d

cd-hit-est: $(COMMON_OBJS) cdhit-est.o
	$(CC) $(CCFLAGS) cdhit-est.o $(COMMON_OBJS) $(LDFLAGS) cd-hit-est

cd-hit-est-2d: $(COMMON_OBJS) cdhit-est-2d.o
	$(CC) $(CCFLAGS) cdhit-est-2d.o $(COMMON_OBJS) $(LDFLAGS) cd-hit-est-2d

cd-hit-div: $(COMMON_OBJS) cdhit-div.o
	$(CC) $(CCFLAGS) cdhit-div.o $(COMMON_OBJS) $(LDFLAGS) cd-hit-div

cd-hit-454: $(COMMON_OBJS) cdhit-454.o
	$(CC) $(CCFLAGS) cdhit-454.o $(COMMON_OBJS) $(LDFLAGS) cd-hit-454

# objects
cdhit-common.o: cdhit-common.c++ cdhit-common.h
	$(CC) $(CCFLAGS) cdhit-common.c++ -c

cdhit-utility.o: cdhit-utility.c++ cdhit-utility.h
	$(CC) $(CCFLAGS) cdhit-utility.c++ -c

cdhit.o: cdhit.c++ cdhit-utility.h
	$(CC) $(CCFLAGS) cdhit.c++ -c

cdhit-2d.o: cdhit-2d.c++ cdhit-utility.h
	$(CC) $(CCFLAGS) cdhit-2d.c++ -c

cdhit-est.o: cdhit-est.c++ cdhit-utility.h
	$(CC) $(CCFLAGS) cdhit-est.c++ -c

cdhit-est-2d.o: cdhit-est-2d.c++ cdhit-utility.h
	$(CC) $(CCFLAGS) cdhit-est-2d.c++ -c

cdhit-div.o: cdhit-div.c++ cdhit-common.h
	$(CC) $(CCFLAGS) cdhit-div.c++ -c

cdhit-454.o: cdhit-454.c++ cdhit-common.h
	$(CC) $(CCFLAGS) cdhit-454.c++ -c

PREFIX ?= /usr/local/bin

install:
	for prog in $(PROGS); do \
		install -m 0755 $$prog $(PREFIX); \
	done
	install -m 0755 *.pl $(PREFIX);
