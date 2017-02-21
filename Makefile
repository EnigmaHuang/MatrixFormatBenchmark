######## DEFINE COMPILER ######################################
# Use the Intel c and C++ compiler
#CC       = icc
#CPP      = icpc

# Use the GNU C and C++ compiler
#CC       = gcc
#CPP      = g++

# Use clang (LLVM) compiler
#CC       = clang
#CPP      = clang++

# PGI
CC		   = pgcc
CPP		   = pgc++

######## DEFINE COMPILER FLAGS ################################
CFLAGS       = -O3
VERBOSEFLAGS = -DVERBOSE

# If you want to make use of the likwid profiler
# uncommend the following two lines
# the evvirement variables LIKWID_LIB and LIKWID_INC has to be set
# on the rrze cluster this variables are already set
LIKWIDFLAGS  = -DUSE_LIKWID $(LIKWID_INC) -DLIKWID_PERFMON
LIKWIDLD_FLAGS = $(LIKWID_LIB) -llikwid -lm

ifeq "$(CC)" "gcc"
	VERBOSEFLAGS += -g -Wall -ansi
	ARCHFLAGS    += -march=native
else ifeq "$(CC)" "icc"
	VERBOSEFLAGS += -g -Wall -ansi
	ARCHFLAGS    += -xhost
else ifeq "$(CC)" "clang"
	VERBOSEFLAGS += -g -Wall -ansi
else ifeq "$(CC)" "pgcc"
#	ARCHFLAGS    += -tp=sandybridge
	VERBOSEFLAGS += -gopt -Minfo=accel,loop,opt,unified,vect,lre,par
endif

CFLAGS    += $(ARCHFLAGS)
CFLAGS    += $(VERBOSEFLAGS)

omp: CFLAGS += -fopenmp $(LIKWIDFLAGS)
acc: CFLAGS += -acc -ta=tesla
#TODO managed, GPU genauer angeben

CPPFLAGS   = $(CFLAGS) -std=c++11
LDFLAGS    = $(LIKWIDLD_FLAGS)
RM         = rm -f

######## DEFINE DEPENDANCY AND RULES ##########################
BIN              = benchmark
OFILES_test      = mmio/mmio.o MMreader.o CSRMatrix.o timing/timing.o test.o
OFILES_benchmark = mmio/mmio.o MMreader.o CSRMatrix.o timing/timing.o benchmark.o

.PHONY: all clean

all: omp

omp: $(BIN)

acc: $(BIN)

clean:
	$(RM) $(BIN) $(OFILES_test) $(OFILES_benchmark)


########## BIN ################################################
test:     $(OFILES_test)
	$(CPP) $(CPPFLAGS) -o $@ $^ $(LDFLAGS)

benchmark:     $(OFILES_benchmark)
	$(CPP) $(CPPFLAGS) -o $@ $^ $(LDFLAGS)

########## C FILES ###########################################
.c.o:
	$(CC) $(CFLAGS) -c $<

mmio/mmio.o: mmio/mmio.c mmio/mmio.h
	$(CC) $(CFLAGS) -c -o $@ $<
timing/timing.o: timing/timing.c timing/timing.h
	$(CC) $(CFLAGS) -c -o $@ $<

########## CPP FILES #########################################
.cpp.o:
	$(CPP) $(CPPFLAGS) -c $<

########## DEPENDENCIES ######################################
test.o: test.cpp CSRMatrix.hpp SellCSigma.hpp MMreader.hpp spMV.hpp
MMreader.o: MMreader.cpp MMreader.hpp mmio/mmio.h
CSRMatrix.o: CSRMatrix.cpp CSRMatrix.hpp MMreader.hpp
benchmark.o: benchmark.cpp CSRMatrix.hpp SellCSigma.hpp MMreader.hpp spMV.hpp
