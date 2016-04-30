#include "CSRMatrix.hpp"

#include <tuple>
#include <iostream>

extern "C"
{
#include <likwid.h>
#include "timing/timing.h"
}

/*****Implementation*CSR_MATRIX***********************************************/
CSR_Matrix::CSR_Matrix( MMreader mmMatrix )
:M_(mmMatrix.getRows()), N_(mmMatrix.getCols()), nz_(mmMatrix.getNonZeros())
//,colInd_(nz_), rowPtr_(M_+1), val_(nz_)
{
    // allocate memory
    colInd_.reserve(nz_);
    val_.reserve(nz_);
    rowPtr_.reserve(M_+1);

    // sort input Matrix by row to create CSR Format
    if( !mmMatrix.isRowSorted() )
        sortByRow(mmMatrix);

    // convert input Format to csr format
    std::vector< std::tuple<int,int,double> > & mmData = mmMatrix.getMatrx();
    int currentRow = -1;
    int currentID  = 0;
    for (auto it=mmData.begin(); it!=mmData.end(); ++it, ++currentID)
    {
        int row = std::get<0>(*it);
        int col = std::get<1>(*it);
        double val = std::get<2>(*it);

        val_.push_back(val);
        colInd_.push_back(col);

        if (row > currentRow)
        {
            rowPtr_.push_back(currentID);
            currentRow = row;
        }
    }
    rowPtr_.push_back(currentID);

/*
    std::cout << "CSR constructed:\n"
              << "(" << M_ << "," << N_ << ") " << nz_ << ":  "
              << colInd_.size() << " " << val_.size() << " "
              << rowPtr_.size() << std::endl;
*/
}

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

#pragma omp for
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
