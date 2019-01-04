DEPENDS := .depend.mk
CXXFLAGS := -g -Wall -O0
SOURCES := count_kmer.cpp
CXX := $(HOME)/local/bin/g++

.PHONY: depend clean

all: count_kmer


count_kmer: count_kmer.o
	$(CXX) $(LDFLAGS) -o $@ $^

depend:
	$(CXX) -MM $(SOURCES) > $(DEPENDS)

clean:
	-rm *.o

-INCLUDE: $(DEPENDS)

