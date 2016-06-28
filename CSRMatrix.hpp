#ifndef CSRMATRIX_HPP
#define CSRMATRIX_HPP

#include <vector>

#include "MMreader.hpp"


/*****Class CSR_MATRIX********************************************************/
class CSR_Matrix
{
public:
    CSR_Matrix( MMreader mmMatrix );        // constructor
    ~CSR_Matrix();                          // destructor

    int            getRows() const     { return M_; }
    int            getCols() const     { return N_; }
    int            getNonZeros() const { return nz_; }
    int const *    getColInd() const   { return colInd_; }
    int const *    getRowPtr() const   { return rowPtr_; }
    double const * getValues() const   { return val_; }

    // We do not need copy and move symantic for this benchmark
    CSR_Matrix(CSR_Matrix const & other) = delete;   // copy constructor
    CSR_Matrix(CSR_Matrix && other) = delete;        // move constructor
    CSR_Matrix & operator= (CSR_Matrix const & other) = delete;  // copy assignment
    CSR_Matrix & operator= (CSR_Matrix && other) = delete;       // move assignment

private:
    int M_, N_, nz_;        // numer of rows, columns and non zeros
    int *colInd_, *rowPtr_; // colum Indices of matrix elements, row Pointer
    double *val_;           // values of matrix elements
    // NOTE: We use row pointer here to ensure NUMA awareness.
};


/*****Free Functions*CSR_MATRIX***********************************************/
/**
 * output operator
 * prints the CSR matrix to os
 */
std::ostream& operator<<(std::ostream& os, CSR_Matrix const & matrix);

#endif
