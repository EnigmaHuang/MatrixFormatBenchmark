#ifndef CRSMATRIX_HPP
#define CSRMatrix_HPP

#include <vector>
#include <tuple>

#include "MMreader.hpp"



/*****Class CSR_MATRIX********************************************************/
class CSR_Matrix
{
public:
    CSR_Matrix() = delete ;
    explicit CSR_Matrix( MMreader mmMatrix );

    int getRows() const { return M_; }
    int getCols() const { return N_; }
    int getNonZeros() const { return nz_; }
    int const * getColInd() const  { return colInd_.data(); }
    int const * getRowPtr() const  { return rowPtr_.data(); }
    double const * getValues() const  { return val_.data(); }

private:
    int M_, N_, nz_;
    std::vector<int> colInd_, rowPtr_;
    std::vector<double> val_;
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
