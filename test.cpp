#include "MMreader.hpp"
#include "CSRMatrix.hpp"
#include "SellCSigma.hpp"
#include "spMV.hpp"

#include <iostream>
#include <vector>
#include <tuple>
#include <cassert>

int main(int argc, char *argv[])
{

/****TEST: IDEBTITY MAXTRIX***************************************************/
    {
    std::cout << "********Identity TEST *****************************" << std::endl;
    // read smal identity
    //CSR_Matrix identity_csr ("matrices/matrix_identity_klein.mtx");
    MMreader identity ("matrices/matrix_identity_klein.mtx");
    CSR_Matrix identity_csr (identity);

    std::vector<double> y,x;
    
    for (int i=0; i<identity_csr.getRows(); ++i)
    {
        x.push_back(i+1);
        y.push_back(42);
    }

     spMV( identity_csr, x.data(), y.data() );

    assert (x == y);
    //std::cout << x;

    std::cout << "Identity_CSR: sucses!" << std::endl;




    SellCSigma_Matrix<3> identity_sell_1_1 (identity,4);

    std::vector<double> m,n;
    
    for (int i=0; i<identity_sell_1_1.getRows(); ++i)
    {
        m.push_back(i+1);
        n.push_back(42);
    }

    spMV( identity_sell_1_1, m.data(), n.data() );

    //std::cout << m;
    //std::cout << n;
    assert (m == n);

    std::cout << "Identity_sell-3-4: sucses!" << std::endl;
    }

/****TEST: MORE COMPLEX MAXTRIX***********************************************/
    {
    std::cout << "********Band TEST *****************************" << std::endl;
    MMreader band ("matrices/matrix_band_klein.mtx");
    CSR_Matrix band_csr (band);

    std::vector<double> y,x;
    
    for (int i=0; i<band_csr.getRows(); ++i)
    {
        x.push_back(1);
        y.push_back(42);
    }

    spMV( band_csr, x.data(), y.data() );

    //assert (y == x);
    //assert (y == 1.);
    //std::cout << x;

    std::cout << "Band CSR: sucses!" << std::endl;




    SellCSigma_Matrix<2> band_sell (band,4);

    std::vector<double> m,n;
    
    for (int i=0; i<band_sell.getRows(); ++i)
    {
        m.push_back(1);
        n.push_back(42);
    }

    spMV( band_sell, m.data(), n.data() );

    assert (m == n);
    assert (n == 1.);
    //std::cout << x;

    std::cout << "Band Sell-2-4: sucses!" << std::endl;
    }


/****TEST: SYMMETRIC MAXTRIX*********************************************/
    {
    std::cout << "********Symetric TEST *****************************" << std::endl;
    MMreader bandSym ("matrices/matrix_band_symmetric_klein.mtx");
    CSR_Matrix band_sym_csr (bandSym);

    std::vector<double> y,x;
    
    for (int i=0; i<band_sym_csr.getRows(); ++i)
    {
        x.push_back(1);
        y.push_back(42);
    }

    spMV( band_sym_csr, x.data(), y.data() );

    assert (y == x);
    assert (y == 1.);
    //std::cout << x;

    std::cout << "Symmetric CSR: sucses!" << std::endl;
    }


/****TEST: This Test must fail!***********************************************/
    {
    std::cout << "********Broken Matrix TEST *****************************" << std::endl;
    MMreader brockenBand ("matrices/matrix_brockenBand_klein.mtx");
    CSR_Matrix brockenBand_csr (brockenBand);

    std::vector<double> y,x;
    
    for (int i=0; i<brockenBand_csr.getRows(); ++i)
    {
        x.push_back(1);
        y.push_back(42);
    }

    spMV( brockenBand_csr, x.data(), y.data() );

    if ( !(y==x) && !(y==1.) )
    {
        std::cout << "brockenBand CSR: sucses!" << std::endl;
    }




    SellCSigma_Matrix<4> brockenBand_sell (brockenBand, 100);

    std::vector<double> m,n;
    
    for (int i=0; i<brockenBand_csr.getRows(); ++i)
    {
        m.push_back(1);
        n.push_back(42);
    }

    spMV( brockenBand_sell, m.data(), n.data() );

    if ( !(m==n) && !(n==1.) )
    {
        std::cout << "brockenBand Sell-4-100: sucses!" << std::endl;
    }


    }

    return 0;
}

