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
 * y=A*x
 * using the CSR Format
 * y and x musst be allocated and valid
 */
template<bool PLUSy=false>
void spMV( CSR_Matrix const & A,
           double const *x,
           double *y
         )
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
    #pragma omp parallel for schedule(runtime)
#endif
    for (int rowID=0; rowID<numRows; ++rowID)
    {
        double tmp = 0.;

        // loop over all elements in row
        for (int rowEntry=rowPtr[rowID]; rowEntry<rowPtr[rowID+1]; ++rowEntry)
        {
            tmp += val[rowEntry] * x[ colInd[rowEntry] ];
        }

        y[rowID] = tmp;
    }

#ifdef USE_LIKWID
    LIKWID_MARKER_STOP("SpMV_CSR");
#endif
}



/*****SELL-C-SIGMA************************************************************/

/**
 * sparse Matrix-Vector multiplication
 * y=A*x
 * using the Sell-C-Sigma Format
 * y and x musst be allocated and valid
 *
 * y musst be large enough to hold values for all paddded rows!
 * x must be permutaed!
 * y will be permutaed!
 */
template< int C>
void spMV( SellCSigma_Matrix const & A,
           double const * x,
           double * y
         )
{
    double const * val       = A.getValues();
    int const * chunkPtr     = A.getChankPtr();
    int const * chunkLength  = A.getChankLength();
    int const * colInd       = A.getColInd();
    int const numRows        = A.getRows();
    int const nonZeros       = A.getNonZeros();
    int const numberOfChunks = A.getNumberOfChunks();
    int const chunkSize      = C;

#ifdef USE_LIKWID
    LIKWID_MARKER_THREADINIT;
    LIKWID_MARKER_START("SpMV_Sell-C-sigma");
#endif

#ifdef _OPENMP
    #pragma omp parallel for schedule(runtime)
#endif
    // loop over all chunks
    for (int chunk=0; chunk < numberOfChunks; ++chunk)
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
                 //cRow<chunkSize && rowID<numRows;
                 cRow<chunkSize;
               ++cRow,           ++rowID
            )
        {
            y[rowID] = tmp[cRow];
        }
    }

#ifdef USE_LIKWID
        LIKWID_MARKER_STOP("SpMV_Sell-C-sigma");
#endif
}

/*wrapper function for dynamic dispatshing*/
void spMV( SellCSigma_Matrix const & A,
           double const * x,
           double * y
         )
{
    int C = A.getChunkSize();

    if (1 == C)
        return spMV<1>(A,x,y);
    else if (2 == C)
        return spMV<2>(A,x,y);
    else if (4 == C)
        return spMV<4>(A,x,y);
    else if (16 == C)
        return spMV<16>(A,x,y);
    else if (32 == C)
        return spMV<32>(A,x,y);
#ifdef SET_C
    else if (SET_C == C)
        return spMV<SET_C>(A,x,y);
#endif
    else
    {
        std::cout << "spMV Kernel for C="<< C << " is not compiled."
                  << " Use 'SET_C=C' as compile time flag to creat this function."
                  << "\nC=1 is used as a fall back function."
                  << std::endl;
        return spMV<1>(A,x,y);
    }

}


#endif
