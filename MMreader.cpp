#include "MMreader.hpp"

#include <cstdio>
#include <cstdlib>
#include <numeric>
#include <algorithm>

extern "C"
{
#include "mmio/mmio.h"
}


/*****Implementation*MMreader*************************************************/
MMreader::MMreader(char const *fileName)
: isRowSorted_(false), isColSorted_(false)
{
    FILE *f;
    MM_typecode matcode;

    if ( NULL == (f=fopen(fileName, "r")) )
    {
        std::cerr << "Can not open file " << fileName << std::endl;
        exit(1);
    }

    // read header of matrix file
    if ( 0 != mm_read_banner(f, &matcode) )
    {
        printf("Could not process Matrix Market banner.\n");
        exit(1);
    }

    // Test propertis 
    if (   !mm_is_matrix(matcode)
        || !mm_is_sparse(matcode)
        || !mm_is_coordinate(matcode)
       )
    {
        printf("Sorry, this application does not support ");
        printf("Market Market type: [%s]\n", mm_typecode_to_str(matcode));
        exit(1);
    }
    if ( mm_is_complex(matcode))
    {
        printf("Complex matrices are not fully suported.");
        printf("The imaginary part is ignored!");
    }

    //TODO use uint

    // get matrix size
    if ( 0 != mm_read_mtx_crd_size(f, &M_, &N_, &nz_) )
    {
        std::cerr << "error while reding size of matrix" << std::endl;
        exit(1);
    }

    // detect symetry and allocate memmory
    if ( mm_is_general(matcode) )
    {
        isSymmetric_ = false;
        matrix_.reserve(nz_);
    }
    else if ( mm_is_symmetric(matcode)
            ||mm_is_hermitian(matcode)
            ||mm_is_skew(matcode)
            )
    {
        isSymmetric_ = true;
        // note: number of nonzeros <= nz_ * 2 (diagonale)
        matrix_.reserve(nz_*2);
    }
    else
    {
        printf("Sorry, this application does not support ");
        printf("Market Market type: [%s]\n", mm_typecode_to_str(matcode));
        exit(1);
    }

    // read matrix
    int row, col;
    double val;
    for (int i=0; i<nz_; ++i)
    {
        if (mm_is_pattern(matcode))
        {
            fscanf(f, "%d %d\n", &row, &col);
            val = 1.;
        }
        else if ( mm_is_complex(matcode))
        {
            fscanf(f, "%d %d %lg %*g\n", &row, &col, &val);
            //TODO fully suport
        }
        else
            fscanf(f, "%d %d %lg\n", &row, &col, &val);

        // adjust from one-baed to zero-based indes
        --row;
        --col;

        matrix_.emplace_back( std::forward_as_tuple(row, col, val) );

        if (isSymmetric_ && (row!=col))
        {
            matrix_.emplace_back( std::forward_as_tuple(col, row, val) );
        }
    }

/*
    std::cout << "MMreader constructed:\n"
                << "(" << M_ << "," << N_ << ") " << nz_ << ": "
                << matrix_.size() << " "
                << isRowSorted_ << " " << isColSorted_ 
                << std::endl;
*/
}


/*****Free Functions MMreader*************************************************/
void sortByRow(MMreader& mmMatrix)
{
    std::vector< std::tuple<int,int,double> > & matrix = mmMatrix.getMatrx(); 

    // first sort by coll
    std::sort( matrix.begin(), matrix.end(),
                [](std::tuple<int,int,double> const &a,
                    std::tuple<int,int,double> const &b)
                    {return std::get<1>(a) < std::get<1>(b);}
                );
    // then (stable!) sort by row
    std::stable_sort( matrix.begin(), matrix.end(),
                        [](std::tuple<int,int,double> const &a,
                            std::tuple<int,int,double> const &b)
                        {return std::get<0>(a) < std::get<0>(b);}
                    );

    mmMatrix.isRowSorted(true);
}

std::vector<int> getValsPerRow(MMreader const & mmMatrix)
{
    std::vector<int> valsPerRow(mmMatrix.getRows(), 0);

    std::vector< std::tuple<int,int,double> > const & mmData = mmMatrix.getMatrx();

    for (auto it=mmData.begin(); it!=mmData.end(); ++it)
    {
        int row = std::get<0>(*it);

        ++valsPerRow[row];
    }

    return valsPerRow;
}

std::vector< std::tuple<int,int> > getRowLengths(MMreader const & mmMatrix)
{
    std::vector< std::tuple<int,int> > rowLengths;
    rowLengths.reserve(mmMatrix.getRows());

    std::vector<int> valuesPerRow = getValsPerRow(mmMatrix);

    for (int row=0; row<mmMatrix.getRows(); ++row)
    {
        rowLengths.emplace_back( std::forward_as_tuple(row, valuesPerRow[row]) );
    }

    return rowLengths;
}

std::vector<int> getOffsets(std::vector<int> const & valuesPerRow)
{
    std::vector<int> offsets;
    offsets.reserve(valuesPerRow.size() + 1);

    offsets.push_back(0);

    std::partial_sum(valuesPerRow.begin(), valuesPerRow.end(),
                     std::back_inserter(offsets));

    return offsets;
}





/*****Output Funktions**********************************/
std::ostream& operator<<( std::ostream& os, std::tuple<int,int> data )
{
    os << "(" << std::get<0>(data) << "," << std::get<1>(data) << ") ";

    return os;
}

std::ostream& operator<<( std::ostream& os, std::tuple<int,int,double> data )
{
    os << "(" << std::get<0>(data) << "," << std::get<1>(data) << ") "
       << std::get<2>(data);

    return os;
}
