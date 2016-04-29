CPP      = g++
CFLAGS   = -O3 -Wall -Winline -Wshadow -ansi -std=c++11
RM       = rm -f

BIN            = example test
OFILES_example = mmio/mmio.o example_read.o
OFILES_test    = mmio/mmio.o MMreader.o CSRMatrix.o test.o


.PHONY: all clean

all: $(BIN)

clean:
	$(RM) $(BIN) $(OFILES_example) $(OFILES_test)

example: $(OFILES_example)
	$(CPP) $(CFLAGS) -o $@ $^

test: $(OFILES_test)
	$(CPP) $(CFLAGS) -o $@ $^

.cpp.o:
	$(CPP) $(CFLAGS) -c $<

mmio/mmio.o: mmio/mmio.c mmio/mmio.h
example_read.o: example_read.c mmio/mmio.h
test.o: test.cpp mmio/mmio.h CSRMatrix.hpp MMreader.hpp
MMreader.o: MMreader.cpp mmio/mmio.h
CSRMatrix.o: CSRMatrix.cpp MMreader.hpp
