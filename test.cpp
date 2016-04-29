#include "MMreader.hpp"
#include "CSRMatrix.hpp"

#include <iostream>
#include <vector>
#include <cassert>

int main(int argc, char *argv[])
{
    //if (argc < 2)
    //{
        //fprintf(stderr, "Usage: %s [martix-market-filename]\n", argv[0]);
        //exit(1);
    //}

    // read MM matraix from file and create CSR matrix
    //MMreader test (argv[1]);
    //CSR_Matrix csr_test(test);

    //std::cout << "Matrix:"<< std::endl;
    //std::cout << csr_test;


/****TEST: IDEBTITY MAXTRIX***************************************************/
    // read smal identity
    MMreader identity ("matrices/matrix_identity_klein.mtx");
    CSR_Matrix identity_csr (identity);

    std::vector<double> y,x;
    
    for (int i=0; i<identity_csr.getRows(); ++i)
    {
        x.push_back(i+1);
        y.push_back(42);
    }

    double runtime = spMV( identity_csr, x.data(), y.data() );

    assert (x == y);
    //std::cout << x;

    std::cout << "Identity: sucses!" << std::endl;

/****TEST: MORE COMPLEX MAXTRIX***********************************************/

    return 0;
}

