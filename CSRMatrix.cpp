#include "CSRMatrix.hpp"

#include <tuple>
#include <iostream>
#include <iterator>
#include <algorithm>


/*****Implementation*CSR_MATRIX***********************************************/
/**
 * constructor
 */
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
#ifdef _OPENMP
    #pragma omp parallel for schedule(runtime)
#endif
    for (int rowID=0; rowID<getRows(); ++rowID)
    {
        rowPtr_[rowID] = offsets[rowID];

        //loop over all elements in Row
        for (int id=offsets[rowID]; id<offsets[rowID+1]; ++id)
        {
            int    col = std::get<1>(mmData[id]);
            double val = std::get<2>(mmData[id]);
        
            val_[id]    = val;
            colInd_[id] = col;
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


/**
 * destructor
 */
CSR_Matrix::~CSR_Matrix()
{
    delete[] val_;
    delete[] rowPtr_;
    delete[] colInd_;
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
