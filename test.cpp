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
    {
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
    std::cout << "Runtime: " << runtime << std::endl;
    }

/****TEST: MORE COMPLEX MAXTRIX***********************************************/
    {
    MMreader band ("matrices/matrix_band_klein.mtx");
    CSR_Matrix band_csr (band);

    std::vector<double> y,x;
    
    for (int i=0; i<band_csr.getRows(); ++i)
    {
        x.push_back(1);
        y.push_back(42);
    }

    double runtime = spMV( band_csr, x.data(), y.data() );

    assert (y == x);
    assert (y == 1.);
    //std::cout << x;

    std::cout << "Band: sucses!" << std::endl;
    std::cout << "Runtime: " << runtime << std::endl;
    }


/****TEST: This Test must fail!***********************************************/
    {
    MMreader brockenBand ("matrices/matrix_brockenBand_klein.mtx");
    CSR_Matrix brockenBand_csr (brockenBand);

    std::vector<double> y,x;
    
    for (int i=0; i<brockenBand_csr.getRows(); ++i)
    {
        x.push_back(1);
        y.push_back(42);
    }

    double runtime = spMV( brockenBand_csr, x.data(), y.data() );

    if ( !(y==x) && !(y==1.) )
    {
        std::cout << "brockenBand: sucses!" << std::endl;
        std::cout << "Runtime: " << runtime << std::endl;
    }

    }

    return 0;
}

