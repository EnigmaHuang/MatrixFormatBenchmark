#include "MMreader.hpp"
#include "CSRMatrix.hpp"

#include <iostream>
#include <vector>
#include <cassert>

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

    double runtime = spMV( csr_matrix, x.data(), y.data() );

    std::cout << "Runtime: " << runtime << std::endl;

    return 0;
}
