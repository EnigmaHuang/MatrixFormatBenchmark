#ifndef SELLCSIGMA_HPP
#define SELLCSIGMA_HPP

#include <vector>
#include <tuple>
#include <iostream>
#include <iterator>
#include <algorithm>
#include <numeric>
#include <cassert>

#include "MMreader.hpp"



/*****Class SELL-C-Sigam******************************************************/
class SellCSigma_Matrix
{
public:
    SellCSigma_Matrix( MMreader mmMatrix, int C, int sigma );   // constructor
    ~SellCSigma_Matrix();                                       // destructor

    int getChunkSize()           const { return C_; }
    int getSigma()               const { return sigma_; }
    int getRows()                const { return M_; }
    int getPaddedRows()          const { return numberOfChunks_*C_; }
    int getCols()                const { return N_; }
    int getNonZeros()            const { return nz_; }
    int getNumberOfChunks()      const { return numberOfChunks_; }
    int getOverhead()            const { return overhead_; }
    int getCapasety()            const { return capasety_; }
    int const * getColInd()      const { return colInd_; }
    int const * getChankPtr()    const { return chunkPtr_; }
    int const * getChankLength() const { return chunkLength_; }
    int const * getPermutation() const { return permute_; }
    double const * getValues()   const { return val_; }

    // We do not need copy and move symantic for this benchmark
    SellCSigma_Matrix(SellCSigma_Matrix const & other) = delete;   // copy constructor
    SellCSigma_Matrix(SellCSigma_Matrix && other) = delete;        // move constructor
    SellCSigma_Matrix & operator= (SellCSigma_Matrix const & other) = delete;  // copy assignment
    SellCSigma_Matrix & operator= (SellCSigma_Matrix && other) = delete;       // move assignment

private:
    int const C_, sigma_;
    int M_, N_, nz_, numberOfChunks_, overhead_, capasety_;
    int *colInd_, *chunkPtr_, *chunkLength_;
    int *permute_;      // Sell-C-sigma row ID -> orginal row ID
    int *antiPermute_;  // orginal row ID -> Sell row ID
    double* val_;
};


/**
 * constructor
 */
SellCSigma_Matrix::SellCSigma_Matrix( MMreader mmMatrix, int C, int const sigma )
:C_(C), sigma_(sigma)
,M_(mmMatrix.getRows()), N_(mmMatrix.getCols())
,nz_(mmMatrix.getNonZeros()), numberOfChunks_((M_-1)/C+1)
,colInd_(nullptr)
,chunkPtr_(new int[numberOfChunks_]), chunkLength_(new int[numberOfChunks_])
,permute_(new int[M_]), antiPermute_(new int[M_])
,val_(nullptr)
{

    // sort input Matrix by row ID
    if( !mmMatrix.isRowSorted() )
        sortByRow(mmMatrix);


    std::vector< std::tuple<int,int,double> > & mmData = mmMatrix.getMatrx();
    std::vector< std::tuple<int, int> > rowLengths = getRowLengths(mmMatrix);

    // sort sigam chunks by row length
    auto begin = rowLengths.begin();
    auto end   = rowLengths.begin() + getSigma();
    for (; end <= rowLengths.end(); begin += getSigma(), end += getSigma() )
    {
        std::sort(begin, end,
                  [](std::tuple<int,int> const & a, std::tuple<int,int> const & b)
                  {return std::get<1>(a) < std::get<1>(b);}
                 );
    }
    begin -= getSigma();
    std::sort(begin, rowLengths.end(),
                [](std::tuple<int,int> const & a, std::tuple<int,int> const & b)
                {return std::get<1>(a) < std::get<1>(b);}
                );


    // determine chunk length and size
    // and set backword permutation
    std::vector<int> valuesPerChunk( getNumberOfChunks() );
#ifdef _OPENMP
    #pragma omp parallel for schedule(runtime)
#endif
    for (int chunk=0; chunk < getNumberOfChunks(); ++chunk)
    {
        int maxRowLenghth = 0;

        for (int i=0,            row=chunk*getChunkSize();
             i<getChunkSize() && row<getRows();
             ++i,                ++row
            )
        {
            if ( maxRowLenghth < std::get<1>(rowLengths[row]) )
                maxRowLenghth = std::get<1>(rowLengths[row]);

            // set backword permutation
            antiPermute_[ std::get<0>(rowLengths[row]) ] = row;
        }

        chunkLength_[chunk] = maxRowLenghth;
        valuesPerChunk[chunk] = maxRowLenghth * getChunkSize();
    }


    // calculate memory usage and allocate memmory for values and colum IDs
    capasety_ = std::accumulate(std::begin(valuesPerChunk),
                                std::end(valuesPerChunk),
                                0
                               );

    val_    = new double[capasety_];
    colInd_ = new int   [capasety_];

    // calulate memory overhead
    overhead_ = capasety_ - getNonZeros();


    // creat Sell-C-sigma data
    std::vector<int> chunkOffset = getOffsets(valuesPerChunk);
    std::vector<int> rowOffset   = getOffsets(getValsPerRow(mmMatrix));

#ifdef _OPENMP
    #pragma omp parallel for schedule(runtime)
#endif
    for (int chunk=0; chunk < getNumberOfChunks(); ++chunk)
    {
        chunkPtr_[chunk] = chunkOffset[chunk];

        for (int j=0; j<chunkLength_[chunk]; ++j)
        {
            for (int i=0,            row=chunk*getChunkSize();
                 i<getChunkSize();
                 ++i,                ++row
                )
            {
                int    col;
                double val;

                if (row<getRows())
                {
                    // set permutation
                    permute_[row] = std::get<0>(rowLengths[row]);

                    // finde values and collumn index
                    if ( j < std::get<1>(rowLengths[row]) )
                    {   // fill with matrix values
                        int id = rowOffset[ permute_[row] ] + j;

                        val = std::get<2>( mmData[id] );
                        col = std::get<1>( mmData[id] );
                    }
                    else
                    {   // fill chunk with 0
                        val = 0.;
                        col = 0;
                    }
                }
                else
                { // add zero rows to end of matrix fill up last chunk
                    val = 0.;
                    col = 0;
                }

                val_   [chunkPtr_[chunk] + i + j*getChunkSize()] = val;
                colInd_[chunkPtr_[chunk] + i + j*getChunkSize()] = antiPermute_[col];
            }
        }
    }

    /*
    std::cout << "Sell-C-sigma constructed:"
              << "\nC: " << getChunkSize() << " sigma: " << getSigma()
              << "\n(" << getRows() << "," << getCols() << ") " << getNonZeros()
              << ":\n";
    for (int i=0; i<valueMemoryUsage; ++i)
    {
        std::cout << getValues()[i] << " (" << getColInd()[i] << ")\n";
    }
    std::cout << std::endl;
    */
}

// destructor
SellCSigma_Matrix::~SellCSigma_Matrix()
{
    delete[] val_;
    delete[] antiPermute_;
    delete[] permute_;
    delete[] chunkLength_;
    delete[] chunkPtr_;
    delete[] colInd_;
}

/*****Free Functions*CSR_MATRIX***********************************************/



#endif
