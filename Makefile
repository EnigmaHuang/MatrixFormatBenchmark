CC       = icc
CPP      = icpc
CFLAGS   = -O3 -xhost -march=native -Wall -Wshadow -ansi -g
CPPFLAGS = $(CFLAGS) -std=c++11
LDFLAGS  = 
OMPFLAG  = -fopenmp
RM       = rm -f

BIN                 = test_seriell test_omp benchmark_seriell benchmark_omp
OFILES_testSer      = mmio/mmio.o MMreader.o CSRMatrix_ser.o timing/timing.o test.o
OFILES_testOMP      = mmio/mmio.o MMreader.o CSRMatrix_omp.o timing/timing.o test.o
OFILES_benchmarkSer = mmio/mmio.o MMreader.o CSRMatrix_ser.o timing/timing.o benchmark.o
OFILES_benchmarkOMP = mmio/mmio.o MMreader.o CSRMatrix_omp.o timing/timing.o benchmark_omp.o


.PHONY: all clean

all: $(BIN)

clean:
	$(RM) $(BIN) $(OFILES_example) $(OFILES_testSer) $(OFILES_testOMP) $(OFILES_benchmarkSer) $(OFILES_benchmarkOMP)


##########BIN#################################################################
test_seriell: $(OFILES_testSer)
	$(CPP) $(CPPFLAGS)            -o $@ $^ $(LDFLAGS)
test_omp:     $(OFILES_testOMP)
	$(CPP) $(CPPFLAGS) $(OMPFLAG) -o $@ $^ $(LDFLAGS)

benchmark_seriell: $(OFILES_benchmarkSer)
	$(CPP) $(CPPFLAGS)            -o $@ $^ $(LDFLAGS)
benchmark_omp:     $(OFILES_benchmarkOMP)
	$(CPP) $(CPPFLAGS) $(OMPFLAG) -o $@ $^ $(LDFLAGS)



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

CSRMatrix_ser.o: CSRMatrix.cpp CSRMatrix.hpp MMreader.hpp
	$(CPP) $(CPPFLAGS) -c -o $@ $<
CSRMatrix_omp.o: CSRMatrix.cpp CSRMatrix.hpp MMreader.hpp
	$(CPP) $(CPPFLAGS) $(OMPFLAG) -c -o $@ $<

benchmark_omp.o: benchmark.cpp CSRMatrix.hpp SellCSigma.hpp MMreader.hpp
	$(CPP) $(CPPFLAGS) $(OMPFLAG) -c -o $@ $<


##########DEPENDENCIES#######################################################
test.o: test.cpp CSRMatrix.hpp SellCSigma.hpp MMreader.hpp spMV.hpp
benchmark.o: benchmark.cpp MMreader.hpp CSRMatrix.hpp SellCSigma.hpp spMV.hpp
MMreader.o: MMreader.cpp MMreader.hpp mmio/mmio.h
