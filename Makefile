GCC      = gcc
CPP      = g++
CFLAGS   = -O3 -Wall -Winline -Wshadow -ansi
CPPFLAGS = $(CFLAGS) -std=c++11
RM       = rm -f

BIN            = example test
OFILES_example = mmio/mmio.o example_read.o
OFILES_test    = mmio/mmio.o MMreader.o CSRMatrix.o timing/timing.o test.o


.PHONY: all clean

all: $(BIN)

clean:
	$(RM) $(BIN) $(OFILES_example) $(OFILES_test)

example: $(OFILES_example)
	$(CPP) $(CPPFLAGS) -o $@ $^

test: $(OFILES_test)
	$(CPP) $(CPPFLAGS) -o $@ $^

# C compiler
mmio/mmio.o: mmio/mmio.c mmio/mmio.h
	$(GCC) $(CFLAGS) -c -o $@ $<
timing/timing.o: timing/timing.c timing/timing.h
	$(GCC) $(CFLAGS) -c -o $@ $<
example_read.o: example_read.c mmio/mmio.h
	$(GCC) $(CFLAGS) -c -o $@ $<

# CPP compiler
.cpp.o:
	$(CPP) $(CPPFLAGS) -c $<

test.o: test.cpp mmio/mmio.h CSRMatrix.hpp MMreader.hpp
MMreader.o: MMreader.cpp mmio/mmio.h
CSRMatrix.o: CSRMatrix.cpp MMreader.hpp
