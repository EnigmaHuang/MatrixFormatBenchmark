#include "MMreader.hpp"
#include "CSRMatrix.hpp"
#include "SellCSigma.hpp"

#include <iostream>
#include <vector>
#include <cassert>



//TODO loop um kernel damit es ein wneig länger dauert
//TODO performace
//TODO andere Formate

int main(int argc, char *argv[])
{
    if (argc < 2)
    {
        fprintf(stderr, "Usage: %s [martix-market-filename]\n", argv[0]);
        exit(1);
    }

    //read MM matraix from file and create CSR matrix
    MMreader mmMatrix (argv[1]);



    /******CSR*******************************************************/
    CSR_Matrix csr_matrix(mmMatrix);
    int const length = csr_matrix.getRows();

    // create vectors (NUMA awareness!)
    std::vector<double> x,y;
    x.reserve(length);
    y.reserve(length);
#pragma omp parallel for schedule(static)
    for (int i=0; i<length; ++i)
    {
        x.push_back(42.);
        y.push_back(0.);
    }

    //for (int i=0; i<1000; ++i)
        spMV( csr_matrix, x.data(), y.data() );
    auto messerment = spMV( csr_matrix, x.data(), y.data() );

    std::cout << "Runtime CSR: " << std::get<0>(messerment) << "sec "
              << "Performance: " << std::get<1>(messerment) << "Flops/sec"
              << std::endl;

    /******SELL*******************************************************/
    SellCSigma_Matrix<4> sell_matrix(mmMatrix, 4);
    int const length_sell = sell_matrix.getRows();

    // create vectors (NUMA awareness!)
    std::vector<double> m,n;
    m.reserve(length_sell);
    n.reserve(length_sell);
#pragma omp parallel for schedule(static)
    for (int i=0; i<length_sell; ++i)
    {
        m.push_back(42.);
        n.push_back(0.);
    }

    //for (int i=0; i<1000; ++i)
        spMV( sell_matrix, m.data(), n.data() );
    auto messerment_sell = spMV( sell_matrix, m.data(), n.data() );

    std::cout << "Runtime SEll: " << std::get<0>(messerment) << "sec "
              << "Performance: " << std::get<1>(messerment) << "Flops/sec"
              << std::endl;
    return 0;
}
