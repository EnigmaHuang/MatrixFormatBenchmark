#include "MMreader.hpp"
//#include "CSRMatrix.hpp"
#include "SellCSigma.hpp"
#include "spMV.hpp"

#include <omp.h>
#include <iostream>
#include <cassert>
#include <cstdlib>
#include <string>

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

    std::cout   << "Matrix Format Benchmark:"
                << "\n\trevisions: " << revisions
                << "\n\tnumber of threads: " << omp_get_num_threads()
                << "\n\tschedular: (" << schedName << ", " << chunkSize << ")"
                << "\n\tmatrix size: " << mmMatrix.getRows() << "x" << mmMatrix.getCols()
                << "\n\tnumber of nonzeros: " << mmMatrix.getNonZeros()
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
    int const length = csr_matrix.getRows();

    double timeing_start, timeing_end, runtime, cpuTime;

    // create vectors (NUMA awareness!)
    double *x = new double[length];
    double *y = new double[length];

    std::cout << "Starting CSR" << std::endl;

    #pragma omp parallel for schedule(runtime)
    for (int i=0; i<length; ++i)
    {
        x[i] = 1.;
        y[i] = 0.;
    }

    timing(&timeing_start, &cpuTime);

    for (int i=0; i<revisions; ++i)
        spMV( csr_matrix, x, y );

    timing(&timeing_end, &cpuTime);
    runtime = timeing_end - timeing_start;

    int flops = csr_matrix.getNonZeros()*2;
    std::cout << "runtime CSR: " << runtime << " sec."
              << " performance: " << static_cast<double>(flops)*revisions / runtime
              << std::endl;

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
    int const length = sell_matrix.getPaddedRows();

    double timeing_start, timeing_end, runtime, cpuTime;

    // create vectors (NUMA awareness!)
    double *x = new double[length];
    double *y = new double[length];

    std::cout << "Starting Sell-" << C << "-" << sigma << std::endl;

    #pragma omp parallel for schedule(runtime)
    for (int i=0; i<length; ++i)
    {
        x[i] = 1.;
        y[i] = 0.;
    }

    timing(&timeing_start, &cpuTime);

    for (int i=0; i<revisions; ++i)
        spMV( sell_matrix, x, y );

    timing(&timeing_end, &cpuTime);
    runtime = timeing_end - timeing_start;

    int flops       = sell_matrix.getNonZeros()*2 ;
    double overhead = static_cast<double>(sell_matrix.getOverhead()) /
                         (sell_matrix.getNonZeros()+sell_matrix.getOverhead());


    std::cout << "runtime Sell-" << C << "-" << sigma << ": " << runtime << " sec."
              << " performance: " << static_cast<double>(flops)*revisions / runtime
              << " flop/s"
              << " overhead: " << overhead*100 << "%"
              << std::endl;

    delete[] x;
    delete[] y;
    }

#ifdef USE_LIKWID
    LIKWID_MARKER_CLOSE;
#endif
    return 0;
}
