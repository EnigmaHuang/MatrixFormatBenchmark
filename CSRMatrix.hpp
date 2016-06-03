#ifndef CRSMATRIX_HPP
#define CSRMATRIX_HPP

#include <vector>
#include <tuple>

#include "MMreader.hpp"

extern "C"
{
#include <likwid.h>
#include "timing/timing.h"
}


/*****Class CSR_MATRIX********************************************************/
class CSR_Matrix
{
public:
    CSR_Matrix( MMreader mmMatrix );        // constructor
    ~CSR_Matrix();                          // destructor

    int getRows() const { return M_; }
    int getCols() const { return N_; }
    int getNonZeros() const { return nz_; }
    int const * getColInd() const  { return colInd_; }
    int const * getRowPtr() const  { return rowPtr_; }
    double const * getValues() const  { return val_; }

    // We do not need copy and move symantic for this benchmark
    CSR_Matrix(CSR_Matrix const & other) = delete;   // copy constructor
    CSR_Matrix(CSR_Matrix && other) = delete;        // move constructor

    CSR_Matrix & operator= (CSR_Matrix const & other) = delete;  // copy assignment
    CSR_Matrix & operator= (CSR_Matrix && other) = delete;       // move assignment

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
    int const rows     = A.getRows();
    int const nonZeros = A.getNonZeros();

//#pragma omp parallel
    //{ // open paralel region
        LIKWID_MARKER_THREADINIT;
        LIKWID_MARKER_START("SpMV_CSR");

        // loop over all rows
#ifdef _OPENMP
        #pragma omp for schedule(runtime)
#endif
        for (int rowID=0; rowID<rows; ++rowID)
        {
            int id = rowPtr[rowID];

            // set y vec to 0
            //y[rowID] = 0;
            double tmp = 0.;

            // loop over all elements in row
            for (; id<rowPtr[rowID+1]; ++id)
            {
                //y[rowID] += val[id] * x[ colInd[id] ];
                tmp += val[id] * x[ colInd[id] ];
            }

            if(PLUSy)
                y[rowID] = alpha * tmp + beta * y[rowID];
            else
                y[rowID] = alpha * tmp;
        }

        LIKWID_MARKER_STOP("SpMV_CSR");

    //} // close paralel region

}

#endif
