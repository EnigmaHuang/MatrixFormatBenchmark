#ifndef SELLCSIGMA_HPP
#define SELLCSIGMA_HPP

#include <vector>
#include <tuple>

#include "MMreader.hpp"



/*****Class SELL-C-Sigam******************************************************/
//TODO sigma als template parameter?
class SellCSigma_Matrix
{
public:
    SellCSigma_Matrix( MMreader mmMatrix, int sigma, int C ); // constructor
    ~SellCSigma_Matrix();                          // destructor

    int getChunkSize() const { return C_; }
    int getSigma() const { return sigma_; }
    int getRows() const { return M_; }
    int getCols() const { return N_; }
    int getNonZeros() const { return nz_; }
    int const * getColInd() const  { return colInd_; }
    int const * getChankPtr() const  { return chunkPtr_; }
    int const * getChankLength() const  { return chunkLength_; }
    int const * getPermutation() const { return permute_; }
    double const * getValues() const  { return val_; }

    // We do not need copy and move symantic for this benchmark
    SellCSigma_Matrix(SellCSigma_Matrix const & other) = delete;   // copy constructor
    SellCSigma_Matrix(SellCSigma_Matrix && other) = delete;        // move constructor

    SellCSigma_Matrix & operator= (SellCSigma_Matrix const & other) = delete;  // copy assignment
    SellCSigma_Matrix & operator= (SellCSigma_Matrix && other) = delete;       // move assignment

private:
    int const sigma_, C_;
    int M_, N_, nz_;
    int *colInd_, *chunkPtr_, *chunkLength_;
    int *permute_;   // Sell-C-sigma row ID -> orginal row ID
    double* val_;
};


/*****Free Functions*CSR_MATRIX***********************************************/
//std::ostream& operator<<(std::ostream& os, SellCSigma_Matrix const & matrix);

/**
 * sparse Matrix-Vector multiplication
 * y=Ax
 * using the CSR Format
 * y and x musst be allocated and valid
 * OMP parallel
 *
 * returns a tuple containg the  runtime and the performance (flops/time)
 * of the kernel
 */
std::tuple<double,double> spMV( SellCSigma_Matrix const & A, double const *x, double *y );

#endif
