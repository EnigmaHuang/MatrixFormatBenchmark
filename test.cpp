#include "MMreader.hpp"
#include "CSRMatrix.hpp"

#include <iostream>
#include <vector>
#include <tuple>
#include <cassert>

int main(int argc, char *argv[])
{

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

    auto messerment = spMV( identity_csr, x.data(), y.data() );

    assert (x == y);
    //std::cout << x;

    std::cout << "Identity: sucses!" << std::endl;
    std::cout << "Runtime: " << std::get<0>(messerment) << "sec "
              << "Perfomance: " << std::get<1>(messerment) << "Flops/sec"
              << std::endl;
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

    auto messerment = spMV( band_csr, x.data(), y.data() );

    assert (y == x);
    assert (y == 1.);
    //std::cout << x;

    std::cout << "Band: sucses!" << std::endl;
    std::cout << "Runtime: " << std::get<0>(messerment) << "sec "
              << "Performance: " << std::get<1>(messerment) << "Flops/sec"
              << std::endl;
    }

/****TEST: MORE SYMMETRIC MAXTRIX*********************************************/
    {
    MMreader bandSym ("matrices/matrix_band_symmetric_klein.mtx");
    CSR_Matrix band_sym_csr (bandSym);

    std::vector<double> y,x;
    
    for (int i=0; i<band_sym_csr.getRows(); ++i)
    {
        x.push_back(1);
        y.push_back(42);
    }

    auto messerment = spMV( band_sym_csr, x.data(), y.data() );

    assert (y == x);
    assert (y == 1.);
    //std::cout << x;

    std::cout << "Symmetric: sucses!" << std::endl;
    std::cout << "Runtime: " << std::get<0>(messerment) << "sec "
              << "Performance: " << std::get<1>(messerment) << "Flops/sec"
              << std::endl;
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

    auto messerment = spMV( brockenBand_csr, x.data(), y.data() );

    if ( !(y==x) && !(y==1.) )
    {
        std::cout << "brockenBand: sucses!" << std::endl;
        std::cout << "Runtime: " << std::get<0>(messerment) << "sec "
                  << "Performance: " << std::get<1>(messerment) << "Flops/sec"
                  << std::endl;
    }

    }

    return 0;
}

