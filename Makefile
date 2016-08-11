# Use the Intel c and C++ compiler
CC       = icc -xhost
CPP      = icpc -xhost

# Use the GNU C and C++ compiler
#CC       = gcc -march=native
#CPP      = g++ -march=native

# Use clang (LLVM) compiler
#CC       = clang
#CPP      = clang++

# If you want to use likwid uncommend this two lines
# and point with LIKDWID_LIB and LIKWID_INC to your likwid instalation
# example: export LIKWID_INC="-I/mnt/opt/likwid-4.0.0_2.11/include"
# on the rrze cluster this variable is already set
LIKWID_FLAGS = -DUSE_LIKWID $(LIKWID_INC) -DLIKWID_PERFMON 
LIKWIDi_LD_FLAGS = $(LIKWID_LIB) -llikwid -lm

CFLAGS   = -O3 -Wall -ansi -g -fopenmp -DVERBOSE $(LIKWID_FLAGS)
CPPFLAGS = $(CFLAGS) -std=c++11
LDFLAGS  = $(LIKWIDi_LD_FLAGS)
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
