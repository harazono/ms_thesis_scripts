DEPENDS := .depend.mk
CXXFLAGS := -g -Wall -O0 -std=c++11
GTESTFLAGS := -lgtest_main -lgtest -lpthread
SOURCES := count_kmer.cpp
TEST_SOURCES := kmer_library_test.cpp
CXX := $(HOME)/local/bin/g++

.PHONY: depend clean test

all: count_kmer


count_kmer: count_kmer.o
	$(CXX) $(LDFLAGS) -o $@ $^

test: build_test
	./kmer_library_test

build_test: kmer_library_test

kmer_library_test: kmer_library_test.o
	$(CXX) $(LDFLAGS) $(GTESTFLAGS) -o $@ $^

depend:
	$(CXX) -MM $(SOURCES) $(TEST_SOURCES) > $(DEPENDS)

clean:
	-rm *.o

-INCLUDE: $(DEPENDS)
