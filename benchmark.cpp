#include "MMreader.hpp"
#include "CSRMatrix.hpp"
#include "SellCSigma.hpp"

#include <iostream>
#include <vector>
#include <cassert>



//TODO performace
//TODO andere Formate
//TODO MEMory alignment

int main(int argc, char *argv[])
{
    if (argc < 2)
    {
        fprintf(stderr, "Usage: %s [martix-market-filename]\n", argv[0]);
        exit(1);
    }

    //read MM matraix from file and create CSR matrix
    MMreader mmMatrix (argv[1]);


    int const revisions = 100000;

    /******CSR*******************************************************/
    {
    CSR_Matrix csr_matrix(mmMatrix);
    int const length = csr_matrix.getRows();

    double timeing_start, timeing_end, runtime, cpuTime;

    // create vectors (NUMA awareness!)
    std::vector<double> x,y;
    x.reserve(length);
    y.reserve(length);
#ifdef _OPENMP
    #pragma omp parallel for schedule(static)
#endif
    for (int i=0; i<length; ++i)
    {
        x.push_back(42.);
        y.push_back(0.);
    }

#ifdef _OPENMP
    #pragma omp parallel
#endif
    {
    timing(&timeing_start, &cpuTime);

    for (int i=0; i<revisions; ++i)
        spMV( csr_matrix, x.data(), y.data() );

    timing(&timeing_end, &cpuTime);
    runtime = timeing_end - timeing_start;
    }

    long flops = csr_matrix.getNonZeros()*2;
    std::cout << "runtime CSR: " << runtime << " sec."
              << " performance: " << flops*revisions / runtime
              << std::endl;

    //auto messerment = spMV( csr_matrix, x.data(), y.data() );

    //std::cout << "Runtime CSR: " << std::get<0>(messerment) << "sec "
              //<< "Performance: " << std::get<1>(messerment) << "Flops/sec"
              //<< std::endl;
    }

    /******SELL*******************************************************/
    {

    SellCSigma_Matrix<4> sell_matrix(mmMatrix, 128);
    int const length_sell = sell_matrix.getRows();

    double timeing_start, timeing_end, runtime, cpuTime;

    // create vectors (NUMA awareness!)
    //std::vector<double> m,n;
    //m.reserve(length_sell);
    //n.reserve(length_sell);
    double *m = (double*) _mm_malloc(sizeof(double)*length_sell,64);
    double *n = (double*) _mm_malloc(sizeof(double)*length_sell,64);
#ifdef _OPENMP
    #pragma omp parallel for schedule(static)
#endif
    for (int i=0; i<length_sell; ++i)
    {
        //m.push_back(42.);
        //n.push_back(0.);
        m[i] = 42.;
        n[i] = 0.;
    }

#ifdef _OPENMP
    #pragma omp parallel
#endif
    {
    timing(&timeing_start, &cpuTime);

    for (int i=0; i<revisions; ++i)
        //spMV( sell_matrix, m.data(), n.data() );
        spMV( sell_matrix, m, n );

    timing(&timeing_end, &cpuTime);
    runtime = timeing_end - timeing_start;
    }

    long long flops    = sell_matrix.getNonZeros()*2 ;
    long long flopsWOh = (sell_matrix.getNonZeros()+sell_matrix.getOverhead())*2 ;

    std::cout << "runtime Sell-4-16: " << runtime << " sec.:"
              << " performance: " << flops*revisions / runtime
              << " performance with overhead: " << flopsWOh*revisions / runtime
              << std::endl;

    //auto messerment_sell = spMV( sell_matrix, m.data(), n.data() );

    //std::cout << "Runtime SEll: " << std::get<0>(messerment) << "sec "
              //<< "Performance: " << std::get<1>(messerment) << "Flops/sec"
              //<< std::endl;
    
    _mm_free(m);
    _mm_free(n);
    }
    return 0;
}
