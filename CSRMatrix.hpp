#ifndef CRSMATRIX_HPP
#define CSRMATRIX_HPP

#include <vector>

#include "MMreader.hpp"

extern "C"
{
//#include <likwid.h>
}


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
 * sparse Matrix-Vector multiplication
 * y=A*x + beta*y
 * using the CSR Format
 * y and x musst be allocated and valid
 * if _OPEMP is set you have to call it inside a OMP parallel region!
 */
template<bool PLUSy=false>
void spMV( CSR_Matrix const & A,
           double const *x,
           double *y,
           double alpha=1.,
           double beta=0.)
{
    double const *val  = A.getValues();
    int const *colInd  = A.getColInd();
    int const *rowPtr  = A.getRowPtr();
    int const numRows  = A.getRows();
    int const nonZeros = A.getNonZeros();

    //LIKWID_MARKER_THREADINIT;
    //LIKWID_MARKER_START("SpMV_CSR");

    // loop over all rows
#ifdef _OPENMP
    #pragma omp for schedule(runtime)
#endif
    for (int rowID=0; rowID<NumRows; ++rowID)
    {
        double tmp = 0.;

        // loop over all elements in row
        for (int rowEntry=rowPtr[rowID]; rowEntry<rowPtr[rowID+1]; ++rowEntry)
        {
            tmp += val[rowEntry] * x[ colInd[rowEntry] ];
        }

        if(PLUSy)
            y[rowID] = alpha * tmp + beta * y[rowID];
        else
        {
            y[rowID] = alpha * tmp;
        }
    }

    //LIKWID_MARKER_STOP("SpMV_CSR");
}


/**
 * output operator
 * prints the CSR matrix to os
 */
std::ostream& operator<<(std::ostream& os, CSR_Matrix const & matrix);

#endif
