#ifndef spMV_HPP
#define spMV_HPP

#include "SellCSigma.hpp"
#include "CSRMatrix.hpp"

#ifdef USE_LIKWID
extern "C"
{
#include <likwid.h>
}
#endif


/*****CSR_MATRIX**************************************************************/
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

#ifdef USE_LIKWID
    LIKWID_MARKER_THREADINIT;
    LIKWID_MARKER_START("SpMV_CSR");
#endif

    // loop over all rows
#ifdef _OPENMP
    #pragma omp for schedule(runtime)
#endif
    for (int rowID=0; rowID<numRows; ++rowID)
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

#ifdef USE_LIKWID
    LIKWID_MARKER_STOP("SpMV_CSR");
#endif
}



/*****SELL-C-SIGMA************************************************************/
/**
 * sparse Matrix-Vector multiplication
 * y=alpha*Ax + beta*y
 * using the CSR Format
 * y and x musst be allocated and valid
 *
 * if _OPEMP is set you have to call it inside a OMP parallel region!
 *
 * x must be permutaed!
 * y will be permutaed!
 */
template< int C, bool PLUSy=false>
void spMV( SellCSigma_Matrix<C> const & A,
           double const * x,
           double * y,
           double alpha=1.,
           double beta=0.)
{
    double const * val       = A.getValues();
    int const * chunkPtr     = A.getChankPtr();
    int const * chunkLength  = A.getChankLength();
    int const * colInd       = A.getColInd();
    int const NumRows        = A.getRows();
    int const nonZeros       = A.getNonZeros();
    int const numberOfChunks = A.getNumberOfChunks();
    int const chunkSize      = C;

#ifdef USE_LIKWID
    LIKWID_MARKER_THREADINIT;
    LIKWID_MARKER_START("SpMV_Sell-C-sigma");
#endif

#ifdef _OPENMP
    #pragma omp for schedule(runtime)
#endif
    // loop over all chunks
    for (int chunk=0; chunk < NumRows/chunkSize; ++chunk)
    {
        int chunkOffset = chunkPtr[chunk];
        double tmp[chunkSize] {};

        // loop over all row elements in chunk
        for (int rowEntry=0; rowEntry<chunkLength[chunk]; ++rowEntry)
        {
            // (auto) vectorised loop over all rows in chunk
            #pragma simd
            for (int cRow=0; cRow<chunkSize; ++cRow)
            {
                tmp[cRow] += val      [chunkOffset + rowEntry*chunkSize + cRow]
                           * x[ colInd[chunkOffset + rowEntry*chunkSize + cRow] ];
            }
        }
        
        // write back result of y = alpha Ax + beta y
        for (int cRow=0,           rowID=chunk*chunkSize;
                 cRow<chunkSize;
               ++cRow,           ++rowID
            )
        {
            if (PLUSy)
                y[rowID] = alpha * tmp[cRow] + beta * y[rowID];
            else
            {
                y[rowID] = alpha * tmp[cRow];
            }
        }
    }

    // loop remainder   -> last (incompleat chunk)
#ifdef _OPENMP
    #pragma omp single
#endif
    if (NumRows/chunkSize != numberOfChunks)
    {
        assert (NumRows/chunkSize == numberOfChunks-1);

        int chunkOffset = chunkPtr[numberOfChunks-1];
        double tmp[chunkSize] {};

        // do MatVecMul
        for (int j=0; j<chunkLength[numberOfChunks-1]; ++j)
        {
            for (int i=0,           row=(numberOfChunks-1)*chunkSize;
                     i<chunkSize && row<NumRows;
                   ++i,           ++row
                )
            {
                tmp[i] += val      [chunkOffset + j*chunkSize + i]
                        * x[ colInd[chunkOffset + j*chunkSize + i] ];
            }
        }
        
        // write back result of y = alpha Ax + beta y
        for (int i=0,           row=(numberOfChunks-1)*chunkSize;
                 i<chunkSize && row<NumRows;
                ++i,          ++row
            )
        {
            if (PLUSy)
                y[row] = alpha * tmp[i] + beta * y[row];
            else
            {
#pragma vector nontemporal
                y[row] = alpha * tmp[i];
            }
        }
    }

#ifdef USE_LIKWID
        LIKWID_MARKER_STOP("SpMV_Sell-C-sigma");
#endif
}


#endif
