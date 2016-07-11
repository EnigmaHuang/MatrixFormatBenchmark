# Use the Intel c and C++ compiler
CC       = icc -xhost
CPP      = icpc -xhost

# Use the GNU C and C++ compiler
#CC       = gcc -march=native
#CPP      = g++ -march=native

# Use clang (LLVM) compiler
#CC       = clang
#CPP      = clang++

CFLAGS   = -O3 -Wall -ansi -g -fopenmp -DVERBOSE -DUSE_LIKWID $(LIKWID_INC) -DLIKWID_PERFMON
CPPFLAGS = $(CFLAGS) -std=c++11
LDFLAGS  =  $(LIKWID_LIB) -llikwid -lm
RM       = rm -f

BIN                 = test_omp benchmark_omp
OFILES_testOMP      = mmio/mmio.o MMreader.o CSRMatrix.o timing/timing.o test.o
OFILES_benchmarkOMP = mmio/mmio.o MMreader.o CSRMatrix.o timing/timing.o benchmark.o


.PHONY: all clean

all: $(BIN)

clean:
	$(RM) $(BIN) $(OFILES_testOMP) $(OFILES_benchmarkOMP)


##########BIN#################################################################
test_omp:     $(OFILES_testOMP)
	$(CPP) $(CPPFLAGS) -o $@ $^ $(LDFLAGS)

benchmark_omp:     $(OFILES_benchmarkOMP)
	$(CPP) $(CPPFLAGS) -o $@ $^ $(LDFLAGS)



##########C FILES############################################################
.c.o:
	$(CC) $(CFLAGS) -c $<

mmio/mmio.o: mmio/mmio.c mmio/mmio.h
	$(CC) $(CFLAGS) -c -o $@ $<
timing/timing.o: timing/timing.c timing/timing.h
	$(CC) $(CFLAGS) -c -o $@ $<

##########CPP FILES##########################################################
.cpp.o:
	$(CPP) $(CPPFLAGS) -c $<


##########DEPENDENCIES#######################################################
test.o: test.cpp CSRMatrix.hpp SellCSigma.hpp MMreader.hpp spMV.hpp
MMreader.o: MMreader.cpp MMreader.hpp mmio/mmio.h
CSRMatrix.o: CSRMatrix.cpp CSRMatrix.hpp MMreader.hpp
benchmark.o: benchmark.cpp CSRMatrix.hpp SellCSigma.hpp MMreader.hpp spMV.hpp
