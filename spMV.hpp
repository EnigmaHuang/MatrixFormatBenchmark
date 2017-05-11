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
 * y= y + A*x
 * using the CSR Format
 * y and x musst be allocated and valid
 * if accelorators are used (openACC) data have to be preent on the device
 */
void spMV( CSR_Matrix const & A,
           double const *x,
           double *y
         )
{
    double const *val     = A.getValues();
    int const *colInd     = A.getColInd();
    int const *rowPtr     = A.getRowPtr();
    int const numRows     = A.getRows();
    int const numNonZeros = A.getNonZeros();

    // loop over all rows
#pragma omp parallel for schedule(runtime)
#pragma acc parallel present(val[0:numNonZeros],            \
                                colInd[0:numNonZeros],      \
                                rowPtr[0:numRows+1],        \
                                x[0:numRows],               \
                                y[0:numRows])               \
                     loop
    for (int rowID=0; rowID<numRows; ++rowID)
    {
        double tmp = y[rowID];

        // loop over all elements in row
        for (int rowEntry=rowPtr[rowID]; rowEntry<rowPtr[rowID+1]; ++rowEntry)
        {
            tmp += val[rowEntry] * x[ colInd[rowEntry] ];
        }

        y[rowID] = tmp;
    }
}



/*****SELL-C-SIGMA************************************************************/

/**
 * sparse Matrix-Vector multiplication
 * y= y + A*x
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
    int const numberOfChunks = A.getNumberOfChunks();
    int const chunkSize      = C;
    int const paddedRows     = A.getPaddedRows();
    int const capacity       = A.getCapasety();

#pragma omp parallel for schedule(runtime)
#pragma acc parallel present(val[0 : capacity],                     \
                             colInd[0 : capacity],                  \
                             chunkPtr[0 : numberOfChunks],          \
                             chunkLength[0 : numberOfChunks],       \
                             x[0 : paddedRows],                     \
                             y[0 : paddedRows])                     \
                     vector_length(C)                               \
            loop
    // loop over all chunks
    for (int chunk=0; chunk < numberOfChunks; ++chunk)
    {
        int chunkOffset = chunkPtr[chunk];
        int rowOffset   = chunk*chunkSize;

        #pragma acc loop vector
        #pragma omp simd
        for (int chunkRow=0; chunkRow<chunkSize; ++chunkRow)
        {
            int globalRow = rowOffset + chunkRow;

            // fill tempory vector with values from y
            double tmp = y[globalRow];

            // loop over all row elements in chunk
            for (int rowEntry=0; rowEntry<chunkLength[chunk]; ++rowEntry)
            {
                tmp +=      val [chunkOffset + rowEntry*chunkSize + chunkRow]
                     * x[ colInd[chunkOffset + rowEntry*chunkSize + chunkRow] ];
            }

            // write back result of y = alpha Ax + beta y
            y[globalRow] = tmp;
        }
    }
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
    else if (64 == C)
      return spMV<64>(A,x,y);
    else if (128 == C)
        return spMV<128>(A,x,y);
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
