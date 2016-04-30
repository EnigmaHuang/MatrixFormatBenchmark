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

    // create vectors
    std::vector<double> x (csr_matrix.getRows(), 42.);
    std::vector<double> y (csr_matrix.getRows(), 0.);

    auto messerment = spMV( csr_matrix, x.data(), y.data() );

    std::cout << "Runtime: " << std::get<0>(messerment) << "sec "
              << "Performance: " << std::get<1>(messerment) << "Flops/sec"
              << std::endl;

    return 0;
}
