#ifndef CRSMATRIX_HPP
#define CSRMatrix_HPP

#include <vector>
#include <tuple>

#include "MMreader.hpp"



/*****Class CSR_MATRIX********************************************************/
class CSR_Matrix
{
public:
    CSR_Matrix( MMreader mmMatrix );        // constructor
    ~CSR_Matrix();                          // destructor

    // We do not need copy and move symantic for this benchmark
    CSR_Matrix(CSR_Matrix const & other) = delete;   // copy constructor
    CSR_Matrix(CSR_Matrix && other) = delete;        // move constructor

    CSR_Matrix & operator= (CSR_Matrix const & other) = delete;  // copy assignment
    CSR_Matrix & operator= (CSR_Matrix && other) = delete;       // move assignment

    int getRows() const { return M_; }
    int getCols() const { return N_; }
    int getNonZeros() const { return nz_; }
    int const * getColInd() const  { return colInd_; }
    int const * getRowPtr() const  { return rowPtr_; }
    double const * getValues() const  { return val_; }

private:
    int M_, N_, nz_;
    int *colInd_, *rowPtr_;
    double* val_;
};


/*****Free Functions*CSR_MATRIX***********************************************/
std::ostream& operator<<(std::ostream& os, CSR_Matrix const & matrix);

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
std::tuple<double,double> spMV( CSR_Matrix const & A, double const *x, double *y );

#endif
