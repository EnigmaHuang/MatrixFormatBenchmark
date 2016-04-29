#include "CSRMatrix.hpp"

#include <tuple>
#include <iostream>


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
std::ostream& operator<<(std::ostream& os, CSR_Matrix const & matrix)
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
