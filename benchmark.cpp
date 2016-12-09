#include "MMreader.hpp"
#include "SellCSigma.hpp"
#include "spMV.hpp"

#include <omp.h>
#include <iostream>
#include <cassert>
#include <cstdlib>
#include <string>
#include <utility>

extern "C"
{
#include "timing/timing.h"

#ifdef USE_LIKWID
#include <likwid.h>
#endif
}


int main(int argc, char *argv[])
{
    if (argc < 2)
    {
        fprintf(stderr, "Usage: %s martix-market-filename [C [sigma [revisions] ] ]\n", argv[0]);
        exit(1);
    }

    //read MM matraix from file and create CSR matrix
    MMreader mmMatrix (argv[1]);


    int revisions = 100;
    if (argc > 4)
        revisions = std::atoi(argv[4]);

#ifdef VERBOSE
#pragma omp parallel
#pragma omp master
{
    int chunkSize;
    omp_sched_t schedType;
    std::string schedName;

    omp_get_schedule(&schedType, &chunkSize);

    switch (schedType)
    {
        case omp_sched_static : schedName = "static"; break;
        case omp_sched_dynamic : schedName = "dynamic"; break;
        case omp_sched_guided : schedName = "guided"; break;
        case omp_sched_auto: schedName = "auto"; break;
    }

    //TODO open ACC infos

    std::cout   << "Matrix Format Benchmark:"
                << "\n\trevisions: " << revisions
                << "\n\tnumber of threads: " << omp_get_num_threads()
                << "\n\tschedular: (" << schedName << ", " << chunkSize << ")"
                << "\n\tmatrix size: " << mmMatrix.getRows() << "x" << mmMatrix.getCols()
                << "\n\tnumber of nonzeros: " << mmMatrix.getNonZeros() << " ("
                    << static_cast<double>(mmMatrix.getNonZeros())/mmMatrix.getRows()
                    << " per Row)"
                << std::endl;
}
#endif


#ifdef USE_LIKWID
    LIKWID_MARKER_INIT;

#pragma omp parallel
{
    LIKWID_MARKER_THREADINIT;
}
#endif

    /******CSR*******************************************************/
    {
    CSR_Matrix csr_matrix(mmMatrix);

    double const *val     = csr_matrix.getValues();
    int const *colInd     = csr_matrix.getColInd();
    int const *rowPtr     = csr_matrix.getRowPtr();
    int const numRows     = csr_matrix.getRows();
    int const length = csr_matrix.getRows();
    int const numNonZeros = csr_matrix.getNonZeros();

    double timeing_start, timeing_end, runtime, cpuTime;

    // create vectors (NUMA awareness!)
    double *x = new double[length];
    double *y = new double[length];

    #pragma omp parallel for schedule(runtime)
    for (int i=0; i<length; ++i)
    {
        x[i] = 1.;
        y[i] = 0.;
    }

    // copy data to device (if necessary)
    #pragma acc data copyin (x[0:length],           \
                             val[0:numNonZeros],    \
                             colInd[0:numNonZeros], \
                             rowPtr[0:numRows+1])   \
                     copyout(y[0:length])
    {
        std::cout << "Starting CSR" << std::endl;

        timing(&timeing_start, &cpuTime);
#ifdef USE_LIKWID
#pragma omp parallel
{
    LIKWID_MARKER_START("SpMV_CSR");
}
#endif

        for (int i=0; i<revisions; ++i)
        {
            spMV( csr_matrix, x, y );
            // swap pointer
            std::swap(x,y);
        }

#ifdef USE_LIKWID
#pragma omp parallel
{
    LIKWID_MARKER_STOP("SpMV_CSR");
}
#endif

        timing(&timeing_end, &cpuTime);
        runtime = timeing_end - timeing_start;

        double flops = numNonZeros * 2;
        std::cout << "runtime CSR: " << runtime << " sec."
                << " performance: " << flops * revisions / runtime
                << std::endl;
    } // copy data back form device

    delete[] x;
    delete[] y;

    }

    /******SELL*******************************************************/
    {

    int C = 4;
    if (argc > 2)
        C = std::atoi(argv[2]);

    int sigma = 1;
    if (argc > 3)
        sigma = std::atoi(argv[3]);


    SellCSigma_Matrix sell_matrix(mmMatrix, C, sigma);

    double const * val       = sell_matrix.getValues();
    int const * colInd       = sell_matrix.getColInd();
    int const * chunkPtr     = sell_matrix.getChankPtr();
    int const * chunkLength  = sell_matrix.getChankLength();
    int const numberOfChunks = sell_matrix.getNumberOfChunks();
    int const length         = sell_matrix.getPaddedRows();
    int const capasety       = sell_matrix.getCapasety();
    int const overhead       = sell_matrix.getOverhead();
    int const nonZeros       = sell_matrix.getNonZeros();

    double timeing_start, timeing_end, runtime, cpuTime;

    // create vectors (NUMA awareness!)
    double *x = new double[length];
    double *y = new double[length];

    #pragma omp parallel for schedule(runtime)
    //TODO like in kernel (junks!!)
    for (int i=0; i<length; ++i)
    {
        x[i] = 1.;
        y[i] = 0.;
    }

    // copy data to device (if necessary)
    #pragma acc data copyin (x[0 : length],                     \
                             val[0 : capasety],                 \
                             colInd[0 : capasety],              \
                             chunkPtr[0 : numberOfChunks],      \
                             chunkLength[0 : numberOfChunks])   \
                     copyout(y[0 : length])
    {

        std::cout << "Starting Sell-" << C << "-" << sigma << std::endl;

        timing(&timeing_start, &cpuTime);

#ifdef USE_LIKWID
#pragma omp parallel
{
        LIKWID_MARKER_START("SpMV_SELL-C-SIGMA");
}
#endif

        for (int i=0; i<revisions; ++i)
        {
            spMV( sell_matrix, x, y );
            // swap pointer
            std::swap(x,y);
        }

#ifdef USE_LIKWID
#pragma omp parallel
{
        LIKWID_MARKER_STOP("SpMV_SELL-C-SIGMA");
}
#endif

        timing(&timeing_end, &cpuTime);
        runtime = timeing_end - timeing_start;

        int flops       = nonZeros * 2;

        std::cout << "runtime Sell-" << C << "-" << sigma << ": " << runtime << " sec."
                  << " performance: " << static_cast<double>(flops)*revisions / runtime
                  << " flop/s overhead: "
                  << static_cast<double>(overhead)/(nonZeros+overhead)*100 << "%"
                  << std::endl;
    } // copy data back from device

    delete[] x;
    delete[] y;
    }

#ifdef USE_LIKWID
    LIKWID_MARKER_CLOSE;
#endif

    return 0;
}
