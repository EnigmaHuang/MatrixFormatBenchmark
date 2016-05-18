#include "SellCSigma.hpp"

#include <tuple>
#include <iostream>
#include <iterator>
#include <algorithm>
#include <numeric>
#include <cassert>

extern "C"
{
#include <likwid.h>
#include "timing/timing.h"
}

/*****Implementation Sell-C-Sigma*********************************************/
SellCSigma_Matrix::SellCSigma_Matrix( MMreader mmMatrix,
                                      int const C, int const sigma )
:sigma_(sigma), C_(C)
,M_(mmMatrix.getRows()), N_(mmMatrix.getCols()), nz_(mmMatrix.getNonZeros())
,colInd_(nullptr)
,chunkPtr_(new int[M_/C]), chunkLength_(new int[M_/C])
,permute_(new int[M_])
,val_(nullptr)
{
    //TODO throw error; no assert
    assert ( 0 == (sigma % getChunkSize()) );
    assert ( 0 == (getRows() % sigma) );

    // sort input Matrix by row ID
    if( !mmMatrix.isRowSorted() )
        sortByRow(mmMatrix);


    std::vector< std::tuple<int,int,double> > & mmData = mmMatrix.getMatrx();
    std::vector< std::tuple<int, int> > rowLengths = getRowLengths(mmMatrix);

    // sort sigam chunks by row length
    for (auto begin = rowLengths.begin(), end = rowLengths.begin()+getSigma();
         end <= rowLengths.end();
         begin += getSigma(), end += getSigma()
        )
    {
        std::sort(begin, end,
                  [](std::tuple<int,int> const & a, std::tuple<int,int> const & b)
                  {return std::get<1>(a) < std::get<1>(b);}
                 );
    }


    // determine chunk length
    std::vector<int> valuesPerChunk( getRows()/getChunkSize() );
#pragma omp parallel for schedule(static)
    for (int chunk=0; chunk < getRows()/getChunkSize(); ++chunk)
    {
        int maxRowLenghth = 0;

        //for (int i=0, row=chun)*C_; i<getChunkSize() && row<getRows(); ++i,++row)
        for (int i=0, row=chunk*getChunkSize(); i<getChunkSize(); ++i,++row)
        {
            if ( maxRowLenghth < std::get<1>(rowLengths[row]) )
                maxRowLenghth = std::get<1>(rowLengths[row]);
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
#pragma omp parallel for schedule(static)
    for (int chunk=0; chunk < getRows()/getChunkSize(); ++chunk)
    {
        chunkPtr_[chunk] = chunkOffset[chunk];

        //TODO unrool (hand oder tamplate(Compiler))
        //TODO loop mit innerer tauschen
        for (int i=0, row=chunk*getChunkSize(); i<getChunkSize(); ++i,++row)
        {
            // set permutation
            permute_[row] = std::get<0>(rowLengths[row]);

            for (int j=0; j<chunkLength_[chunk]; ++j)
            {
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
                colInd_[chunkPtr_[chunk] + i + j*getChunkSize()] = col;
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


SellCSigma_Matrix::~SellCSigma_Matrix()
{
    delete[] val_;
    delete[] permute_;
    delete[] chunkLength_;
    delete[] chunkPtr_;
    delete[] colInd_;
}


/*****Free Functions*CSR_MATRIX***********************************************/
std::tuple<double,double> spMV( SellCSigma_Matrix const & A, double const *x, double *y )
{
    double const *val      = A.getValues();
    int const *colInd      = A.getColInd();
    int const *chunkPtr    = A.getChankPtr();
    int const *chunkLength = A.getChankLength();
    int const *permute     = A.getPermutation();
    int const rows         = A.getRows();
    int const chunkSize    = A.getChunkSize();
    int const nonZeros     = A.getNonZeros();

    double timeing_start, timeing_end, runtime, cpuTime;
    double performance;

#pragma omp parallel
    { // open paralel region
        timing(&timeing_start, &cpuTime);

        LIKWID_MARKER_THREADINIT;
        LIKWID_MARKER_START("SpMV_Sell-C-sigma");

#pragma omp for schedule(static)
    for (int chunk=0; chunk < rows/chunkSize; ++chunk)
    {
        int chunkOffset = chunkPtr[chunk];

        // zero out target vector
        for (int i=0, row=chunk*chunkSize; i<chunkSize; ++i,++row)
        {
            y[ permute[row] ] = 0.;
        }

        // do MatVecMul
        for (int j=0; j<chunkLength[chunk]; ++j)
        {
            //TODO unrool (hand oder tamplate(Compiler))
            for (int i=0, row=chunk*chunkSize; i<chunkSize; ++i,++row)
            {
                y[ permute[row] ] += val      [chunkOffset + j*chunkSize + i]
                                   * x[ colInd[chunkOffset + j*chunkSize + i] ];
            }
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

