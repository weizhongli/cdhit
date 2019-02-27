CXX ?= g++

# without OpenMP

# with OpenMP
# in command line: 
# make openmp=yes
ifeq ($(openmp),no)
  CXXFLAGS += -DNO_OPENMP
else
  CXXFLAGS += -fopenmp
endif

# support debugging
# in command line:
# make debug=yes
# make openmp=yes debug=yes
ifeq ($(debug),yes)
CXXFLAGS += -ggdb
else
CXXFLAGS ?= -O2
endif

ifdef MAX_SEQ
CXXFLAGS += -DMAX_SEQ=$(MAX_SEQ)
endif

#LDFLAGS = -static -lz -o
LDFLAGS += -lz -o

PROGS = cd-hit cd-hit-est cd-hit-2d cd-hit-est-2d cd-hit-div cd-hit-454

# Propagate hardening flags
CXXFLAGS := $(CPPFLAGS) $(CXXFLAGS)

.c++.o:
	$(CXX) $(CXXFLAGS) -c $<

all: $(PROGS)

clean:
	rm -f *.o $(PROGS)

# programs

cd-hit: cdhit-common.o cdhit-utility.o cdhit.o
	$(CXX) $(CXXFLAGS) cdhit.o cdhit-common.o cdhit-utility.o $(LDFLAGS) cd-hit

cd-hit-2d: cdhit-common.o cdhit-utility.o cdhit-2d.o
	$(CXX) $(CXXFLAGS) cdhit-2d.o cdhit-common.o cdhit-utility.o $(LDFLAGS) cd-hit-2d

cd-hit-est: cdhit-common.o cdhit-utility.o cdhit-est.o
	$(CXX) $(CXXFLAGS) cdhit-est.o cdhit-common.o cdhit-utility.o $(LDFLAGS) cd-hit-est

cd-hit-est-2d: cdhit-common.o cdhit-utility.o cdhit-est-2d.o
	$(CXX) $(CXXFLAGS) cdhit-est-2d.o cdhit-common.o cdhit-utility.o $(LDFLAGS) cd-hit-est-2d

cd-hit-div: cdhit-common.o cdhit-utility.o cdhit-div.o
	$(CXX) $(CXXFLAGS) cdhit-div.o cdhit-common.o cdhit-utility.o $(LDFLAGS) cd-hit-div

cd-hit-454: cdhit-common.o cdhit-utility.o cdhit-454.o
	$(CXX) $(CXXFLAGS) cdhit-454.o cdhit-common.o cdhit-utility.o $(LDFLAGS) cd-hit-454

# objects
cdhit-common.o: cdhit-common.c++ cdhit-common.h
	$(CXX) $(CXXFLAGS) cdhit-common.c++ -c

cdhit-utility.o: cdhit-utility.c++ cdhit-utility.h
	$(CXX) $(CXXFLAGS) cdhit-utility.c++ -c

cdhit.o: cdhit.c++ cdhit-utility.h
	$(CXX) $(CXXFLAGS) cdhit.c++ -c

cdhit-2d.o: cdhit-2d.c++ cdhit-utility.h
	$(CXX) $(CXXFLAGS) cdhit-2d.c++ -c

cdhit-est.o: cdhit-est.c++ cdhit-utility.h
	$(CXX) $(CXXFLAGS) cdhit-est.c++ -c

cdhit-est-2d.o: cdhit-est-2d.c++ cdhit-utility.h
	$(CXX) $(CXXFLAGS) cdhit-est-2d.c++ -c

cdhit-div.o: cdhit-div.c++ cdhit-common.h
	$(CXX) $(CXXFLAGS) cdhit-div.c++ -c

cdhit-454.o: cdhit-454.c++ cdhit-common.h
	$(CXX) $(CXXFLAGS) cdhit-454.c++ -c

PREFIX ?= $(DESTDIR)/usr/local/bin

install:
	for prog in $(PROGS); do \
		install -m 0755 $$prog $(PREFIX); \
	done
	install -m 0755 *.pl $(PREFIX);
