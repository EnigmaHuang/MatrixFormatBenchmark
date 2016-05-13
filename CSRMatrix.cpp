#include "CSRMatrix.hpp"

#include <tuple>
#include <iostream>
#include <iterator>
#include <algorithm>

extern "C"
{
#include <likwid.h>
#include "timing/timing.h"
}

/*****Implementation*CSR_MATRIX***********************************************/
CSR_Matrix::CSR_Matrix( MMreader mmMatrix )
:M_(mmMatrix.getRows()), N_(mmMatrix.getCols()), nz_(mmMatrix.getNonZeros())
,colInd_(new int[nz_]), rowPtr_(new int[M_+1]), val_(new double[nz_])
{
    // sort input Matrix by row to create CSR Format
    if( !mmMatrix.isRowSorted() )
        sortByRow(mmMatrix);

    std::vector< std::tuple<int,int,double> > & mmData = mmMatrix.getMatrx();
    std::vector<int> valuesPerRow = getValsPerRow(mmMatrix);
    std::vector<int> offsets      = getOffsets(valuesPerRow);

    // convert input Format to csr format (NUMA awareness!)
#pragma omp parallel for schedule(static)
    for (int rowID=0; rowID<getRows(); ++rowID)
    {
        rowPtr_[rowID] = offsets[rowID];

        //loop over all elements in Row
        for (int i=offsets[rowID]; i<offsets[rowID+1]; ++i)
        {
            int    col = std::get<1>(mmData[i]);
            double val = std::get<2>(mmData[i]);
        
            val_[i]    = val;
            colInd_[i] = col;
        }
    }

    rowPtr_[M_] = offsets[M_];

/*
    std::cout << "CSR constructed:\n"
              << "(" << M_ << "," << N_ << ") " << nz_ << ":  "
              << colInd_.size() << " " << val_.size() << " "
              << rowPtr_.size() << std::endl;
*/
}


CSR_Matrix::~CSR_Matrix()
{
    delete[] colInd_;
    delete[] val_;
    delete[] rowPtr_;
}

/*
CSR_Matrix::CSR_Matrix(CSR_Matrix const & other)
:M_(other.M_), N_(other.N_), nz_(other.nz_),
 colInd_(new int[nz_]), rowPtr_(new int[M_+1]), val_(new double[nz_])
{

    std::copy(other.colInd_, other.colInd_+nz_, colInd_);
    std::copy(other.rowPtr_, other.rowPtr_+M_+1, rowPtr_);
    std::copy(other.val_, other.val_+nz_, val_);
}
*/

/*
CSR_Matrix & CSR_Matrix::operator= (CSR_Matrix const & other)
{
    CSR_Matrix tmp (other);
    std::swap ( M_, tmp.M_ );
    std::swap ( N_, tmp.N_ );
    std::swap ( nz_, tmp.nz_ );
    std::swap ( colInd_, tmp.colInd_ );
    std::swap ( rowPtr_, tmp.rowPtr_ );
    std::swap ( val_, tmp.val_ );

    return *this;
}
*/

/*
CSR_Matrix::CSR_Matrix( CSR_Matrix &&other )
:M_(other.M_), N_(other.N_), nz_(other.nz_),
 colInd_(new int[nz_]), rowPtr_(new int[M_+1]), val_(new double[nz_])
{
    other.M_ = 0;
    other.N_ = 0;
    other.nz_ = 0;
    other.colInd_ = nullptr;
    other.rowPtr_ = nullptr;
    other.val_ = nullptr;
}
*/

/*
CSR_Matrix & CSR_Matrix::operator= (CSR_Matrix && other)
{
    CSR_Matrix tmp (other);
    std::swap ( M_, other.M_ );
    std::swap ( N_, other.N_ );
    std::swap ( nz_, other.nz_ );
    std::swap ( colInd_, other.colInd_ );
    std::swap ( rowPtr_, other.rowPtr_ );
    std::swap ( val_, other.val_ );

    return *this;
}
*/

/*****Free Functions*CSR_MATRIX***********************************************/
std::ostream& operator<<( std::ostream& os, CSR_Matrix const & matrix )
{

    int const * colInd = matrix.getColInd();
    int const * rowPtr = matrix.getRowPtr();
    double const * val = matrix.getValues();

    int row = 0;
    for (int rowID=0; rowID<matrix.getRows(); ++rowID)
    {
        for (int id=rowPtr[rowID]; id<rowPtr[rowID+1]; ++id)
        {
            os  << val[id] << " (" << row
                << ", " << colInd[id] << ")" << std::endl;
        }
        ++row;
    }

    return os;
}

std::tuple<double,double> spMV( CSR_Matrix const & A, double const *x, double *y )
{
    double const *val  = A.getValues();
    int const *colInd  = A.getColInd();
    int const *rowPtr  = A.getRowPtr();
    int const rows     = A.getRows();
    int const nonZeros = A.getNonZeros();

    double timeing_start, timeing_end, runtime, cpuTime;
    double performance;

#pragma omp parallel
    { // open paralel region
        timing(&timeing_start, &cpuTime);

        LIKWID_MARKER_THREADINIT;
        LIKWID_MARKER_START("SpMV_CSR");

#pragma omp for schedule(static)
        // loop over all rows
        for (int rowID=0; rowID<rows; ++rowID)
        {
            int id = rowPtr[rowID];

            // set y vec to 0
            y[rowID] = 0;

            // loop over all elements in row
            for (; id<rowPtr[rowID+1]; ++id)
            {
                y[rowID] += val[id] * x[ colInd[id] ];

            }
        }

        LIKWID_MARKER_STOP("SpMV_CSR");

        timing(&timeing_end, &cpuTime);
        runtime = timeing_end - timeing_start;

        int flops = nonZeros*2 - rows;
        performance = flops/runtime;

    } // close paralel region

    return std::make_tuple(runtime, performance);
}

