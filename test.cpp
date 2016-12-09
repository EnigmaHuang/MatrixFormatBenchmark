#include "MMreader.hpp"
#include "CSRMatrix.hpp"
#include "SellCSigma.hpp"
#include "spMV.hpp"

#include <iostream>
#include <vector>
#include <tuple>
#include <cassert>
#include <numeric>
#include <cmath>

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
        x.push_back(1.);
        y.push_back(42.);
    }

     spMV( identity_csr, x.data(), y.data() );

    //assert (x == y);
    //std::cout << y;
    assert (y == 43. );

    std::cout << "Identity_CSR: sucses!" << std::endl;




    SellCSigma_Matrix identity_sell_1_1 (identity,2,5);

    std::vector<double> m,n;
    
    for (int i=0; i<identity_sell_1_1.getPaddedRows(); ++i)
    {
        m.push_back(1.);
        n.push_back(42.);
    }

    spMV( identity_sell_1_1, m.data(), n.data() );

    //std::cout << m;
    //std::cout << n;
    m.resize(identity_sell_1_1.getRows());
    n.resize(identity_sell_1_1.getRows());
    //assert (m == n);
    assert (n == 43.);

    std::cout << "Identity_sell-2-5: sucses!" << std::endl;
    }

/****TEST: MORE COMPLEX MAXTRIX***********************************************/
    {
    std::cout << "********Band TEST *****************************" << std::endl;
    MMreader band ("matrices/matrix_band_klein.mtx");
    CSR_Matrix band_csr (band);

    std::vector<double> y,x;
    
    for (int i=0; i<band_csr.getRows(); ++i)
    {
        x.push_back(1.);
        y.push_back(0.);
    }

    spMV( band_csr, x.data(), y.data() );

    assert (y == x);
    assert (y == 1.);
    //std::cout << x;

    std::cout << "Band CSR: sucses!" << std::endl;




    SellCSigma_Matrix band_sell (band,2,4);

    std::vector<double> m,n;
    
    for (int i=0; i<band_sell.getPaddedRows(); ++i)
    {
        m.push_back(1.);
        n.push_back(0.);
    }

    spMV( band_sell, m.data(), n.data() );

    m.resize(band_sell.getRows());
    n.resize(band_sell.getRows());
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
        x.push_back(1.);
        y.push_back(0.);
    }

    spMV( band_sym_csr, x.data(), y.data() );

    assert (y == x);
    assert (y == 1.);
    //std::cout << x;

    std::cout << "Symmetric: sucses!" << std::endl;
    }


/****TEST: This Test must fail!***********************************************/
    {
    std::cout << "********Broken Matrix TEST *****************************" << std::endl;
    MMreader brockenBand ("matrices/matrix_brockenBand_klein.mtx");
    CSR_Matrix brockenBand_csr (brockenBand);

    std::vector<double> y,x;
    
    for (int i=0; i<brockenBand_csr.getRows(); ++i)
    {
        x.push_back(1.);
        y.push_back(0.);
    }

    spMV( brockenBand_csr, x.data(), y.data() );

    if ( !(y==x) && !(y==1.) )
    {
        std::cout << "brockenBand CSR: sucses!" << std::endl;
    }
    else
    {
        std::cout << "Faild!" << std::endl;
        return 1;
    }




    SellCSigma_Matrix brockenBand_sell (brockenBand, 4, 100);

    std::vector<double> m,n;
    
    for (int i=0; i<brockenBand_sell.getPaddedRows(); ++i)
    {
        m.push_back(1.);
        n.push_back(0.);
    }

    spMV( brockenBand_sell, m.data(), n.data() );

    m.resize(brockenBand_sell.getRows());
    n.resize(brockenBand_sell.getRows());
    if ( !(m==n) && !(n==1.) )
    {
        std::cout << "brockenBand Sell-4-100: sucses!" << std::endl;
    }
    else
    {
        std::cout << "Faild!" << std::endl;
        return 1;
    }


    }

/****Randome matrix test***********************************************/
    if(argc == 2)
    {
        std::cout << "********Randome Matrix Test*****************************" << std::endl;

        MMreader mmMatrix (argv[1]);
        CSR_Matrix csrMatrix(mmMatrix);
        SellCSigma_Matrix sellMatrix(mmMatrix,4, 128);

        std::vector<double> x,yCSR,ySell;
    
        for (int i=0; i<sellMatrix.getPaddedRows(); ++i)
        {
            x.push_back(1.);
            yCSR.push_back(0.);
            ySell.push_back(0.);
        }

        spMV (csrMatrix, x.data(), yCSR.data());
        spMV (sellMatrix, x.data(), ySell.data());

        yCSR.resize(sellMatrix.getRows());
        ySell.resize(sellMatrix.getRows());

        double csrSum  = std::accumulate(yCSR.begin(),  yCSR.end(),  0.);
        double sellSum = std::accumulate(ySell.begin(), ySell.end(), 0.);

        double epsilon = 1e-6;

        if (std::abs(csrSum - sellSum) <= epsilon)
        {
            std::cout << "Randome matrix test: sucses!" << std::endl;
        }
        else
        {
            std::cout << "Faild!" << std::endl;
            return 1;
        }
    }





    return 0;
}

