#include "MMreader.hpp"
#include "CSRMatrix.hpp"

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

    auto messerment = spMV( csr_matrix, x.data(), y.data() );

    std::cout << "Runtime: " << std::get<0>(messerment) << "sec "
              << "Performance: " << std::get<1>(messerment) << "Flops/sec"
              << std::endl;

    return 0;
}
