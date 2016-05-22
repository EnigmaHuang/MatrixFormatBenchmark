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

extern "C"
{
#include <likwid.h>
#include "timing/timing.h"
}


/*****Class SELL-C-Sigam******************************************************/
//TODO sigma als template parameter?
template <int C>
class SellCSigma_Matrix
{
public:
    SellCSigma_Matrix( MMreader mmMatrix, int sigma ); // constructor
    ~SellCSigma_Matrix();                          // destructor

    int getChunkSize() const { return C; }
    int getSigma() const { return sigma_; }
    int getRows() const { return M_; }
    int getCols() const { return N_; }
    int getNonZeros() const { return nz_; }
    int getNumberOfChunks() const { return numberOfChunks_; }
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
    int const sigma_;
    int M_, N_, nz_, numberOfChunks_;
    int *colInd_, *chunkPtr_, *chunkLength_;
    //TODO neue namen!
    int *permute_;          // Sell-C-sigma row ID -> orginal row ID
    int *permuteMinus1_;    // orginal row ID -> Sell row ID
    double* val_;
};


template <int C>
SellCSigma_Matrix<C>::SellCSigma_Matrix( MMreader mmMatrix,
                                      int const sigma )
:sigma_(sigma)
,M_(mmMatrix.getRows()), N_(mmMatrix.getCols())
,nz_(mmMatrix.getNonZeros()), numberOfChunks_((M_-1)/C+1)
,colInd_(nullptr)
//,chunkPtr_(new int[(M_-1)/C+1]), chunkLength_(new int[(M_-1)/C+1])
,chunkPtr_(new int[numberOfChunks_]), chunkLength_(new int[numberOfChunks_])
,permute_(new int[M_]), permuteMinus1_(new int[M_])
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
    for (;
         end <= rowLengths.end();
         begin += getSigma(), end += getSigma()
        )
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
#pragma omp parallel for schedule(runtime)
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
            permuteMinus1_[ std::get<0>(rowLengths[row]) ] = row;
        }

        chunkLength_[chunk] = maxRowLenghth;
        valuesPerChunk[chunk] = maxRowLenghth * getChunkSize();

    }


    // calculate memory usage and allocate memmory for values and colum IDs
    size_t valueMemoryUsage = std::accumulate(std::begin(valuesPerChunk),
                                              std::end(valuesPerChunk),
                                              0
                                             );
    val_    = new double[valueMemoryUsage];
    colInd_ = new int[valueMemoryUsage];


    std::vector<int> chunkOffset = getOffsets(valuesPerChunk);
    std::vector<int> rowOffset   = getOffsets(getValsPerRow(mmMatrix));

    // creat Sell-C-sigma data
#pragma omp parallel for schedule(runtime)
    for (int chunk=0; chunk < getNumberOfChunks(); ++chunk)
    {
        chunkPtr_[chunk] = chunkOffset[chunk];

        for (int j=0; j<chunkLength_[chunk]; ++j)
        {
            for (int i=0,            row=chunk*getChunkSize();
                 i<getChunkSize() && row<getRows();
                 ++i,                ++row
                )
            {
                // set permutation
                permute_[row] = std::get<0>(rowLengths[row]);

                int    col;
                double val;

                if ( j < std::get<1>(rowLengths[row]) )
                {   // matrix values
                    int id = rowOffset[ permute_[row] ] + j;

                    val = std::get<2>( mmData[id] );
                    col = std::get<1>( mmData[id] );
                }
                else
                {   // fill chunk with 0
                    val = 0.;
                    col = 0; //TODO irgendwas kluges hier? zB row?  (out of bounds?!)
                }

                val_   [chunkPtr_[chunk] + i + j*getChunkSize()] = val;
                colInd_[chunkPtr_[chunk] + i + j*getChunkSize()] = permuteMinus1_[col];
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

template<int C>
SellCSigma_Matrix<C>::~SellCSigma_Matrix()
{
    delete[] val_;
    delete[] permute_;
    delete[] chunkLength_;
    delete[] chunkPtr_;
    delete[] colInd_;
}

/*****Free Functions*CSR_MATRIX***********************************************/
//std::ostream& operator<<(std::ostream& os, SellCSigma_Matrix const & matrix);

/**
 * sparse Matrix-Vector multiplication
 * y=alpha*Ax + beta*y
 * using the CSR Format
 * y and x musst be allocated and valid
 * OMP parallel
 *
 * x must be permutaed!
 * y will be permutaed!
 *
 * returns a tuple containg the  runtime and the performance (flops/time)
 * of the kernel
 */
//TODO rename aAxpby
template< int C, bool PLUSY=false>
std::tuple<double,double> spMV( SellCSigma_Matrix<C> const & A,
                                double const *x,
                                double *y,
                                double alpha=1.,    //TODO reicht das zur runtime wenn das gleich 1 ist?
                                double beta=0.)
{
    double const *val        = A.getValues();
    int const *colInd        = A.getColInd();
    int const *chunkPtr      = A.getChankPtr();
    int const *chunkLength   = A.getChankLength();
    int const rows           = A.getRows();
    int const nonZeros       = A.getNonZeros();
    int const numberOfChunks = A.getNumberOfChunks();
    int const chunkSize      = C;
//#define chunkSize C

    double timeing_start, timeing_end, runtime, cpuTime;
    double performance;

#pragma omp parallel
    { // open paralel region
        timing(&timeing_start, &cpuTime);

        LIKWID_MARKER_THREADINIT;
        LIKWID_MARKER_START("SpMV_Sell-C-sigma");

#pragma omp for schedule(runtime)
    for (int chunk=0; chunk < numberOfChunks; ++chunk)
    {
        int chunkOffset = chunkPtr[chunk];
        double tmp[chunkSize] {};

        // do MatVecMul
        for (int j=0; j<chunkLength[chunk]; ++j)
        {
            //TODO unrool (überprüfen)
            //TODO compiler flag march native?
            for (int i=0,           row=chunk*chunkSize;
                     i<chunkSize && row<rows;
                   ++i,           ++row
                )
            {
                tmp[i] += val      [chunkOffset + j*chunkSize + i]
                        * x[ colInd[chunkOffset + j*chunkSize + i] ];
            }
        }
        
        // write back result of y = alpha Ax + beta y
        for (int i=0,           row=chunk*chunkSize;
                    i<chunkSize && row<rows;
                ++i,           ++row
            )
        {
            if (PLUSY)
                y[row] = alpha * tmp[i] + beta * y[row];
            else
                y[row] = alpha * tmp[i];
        }

    }

        LIKWID_MARKER_STOP("SpMV_Sell-C-sigma");

        timing(&timeing_end, &cpuTime);
        runtime = timeing_end - timeing_start;

        int flops = nonZeros*2 - rows;
        performance = flops/runtime;

    } // close paralel region

    return std::forward_as_tuple(runtime, performance);
}

#endif
