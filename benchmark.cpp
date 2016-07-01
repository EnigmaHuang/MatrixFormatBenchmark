#include "MMreader.hpp"
//#include "CSRMatrix.hpp"
#include "SellCSigma.hpp"
#include "spMV.hpp"

#include <iostream>
#include <cassert>

extern "C"
{
#include "timing/timing.h"
}


int main(int argc, char *argv[])
{
    if (argc < 2)
    {
        fprintf(stderr, "Usage: %s martix-market-filename\n", argv[0]);
        exit(1);
    }

    //read MM matraix from file and create CSR matrix
    MMreader mmMatrix (argv[1]);


    int const revisions = 100;  //TODO comandline parameter

    /******CSR*******************************************************/
    {
    CSR_Matrix csr_matrix(mmMatrix);
    int const length = csr_matrix.getRows();

    double timeing_start, timeing_end, runtime, cpuTime;

    // create vectors (NUMA awareness!)
    double *x = new double[length];
    double *y = new double[length];
#ifdef _OPENMP
    #pragma omp parallel for schedule(runtime)
#endif
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

    //TODO compile for different Cs
    //TODO sigma als comand line parameter
    int C = 4;
    int sigma = 1;

    SellCSigma_Matrix<4> sell_matrix(mmMatrix, sigma);
    int const length = sell_matrix.getPaddedRows();

    double timeing_start, timeing_end, runtime, cpuTime;

    // create vectors (NUMA awareness!)
    double *x = new double[length];
    double *y = new double[length];
#ifdef _OPENMP
    #pragma omp parallel for schedule(runtime)
#endif
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
              << " overhead: " << overhead*100 << "%"
              << std::endl;

    delete[] x;
    delete[] y;
    }

    return 0;
}
