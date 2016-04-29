#include "MMreader.hpp"

#include <cstdio>
#include <cstdlib>

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
        exit(1); //TODO exeptopn -> more cpp style
    }

    // read header of matrix file
    if ( 0 != mm_read_banner(f, &matcode) )
    {
        printf("Could not process Matrix Market banner.\n");
        exit(1);
    }

    // Test propertis 
    if ( !mm_is_matrix(matcode) || !mm_is_sparse(matcode) )
    {
        printf("Sorry, this application does not support ");
        printf("Market Market type: [%s]\n", mm_typecode_to_str(matcode));
        exit(1);
    }
    if ( mm_is_complex(matcode) )
    {
        printf("Complex numbers not suported\n");
        exit(1);
    }

    // get matrix size
    if ( 0 != mm_read_mtx_crd_size(f, &M_, &N_, &nz_) )
    {
        std::cerr << "error while reding size of matrix" << std::endl;
        exit(1);
    }

    // allocate memmory
    matrix_.reserve(nz_);

    // read matrix
    int row, col;
    double val;
    for (int i=0; i<nz_; ++i)
    {
        fscanf(f, "%d %d %lg\n", &row, &col, &val);

        // adjust from one-baed to zero-based indes
        --row;
        --col;

        matrix_.emplace_back( std::forward_as_tuple(row, col, val) );
    }

/*
    std::cout << "MMreader constructed:\n"
                << "(" << M_ << "," << N_ << ") " << nz_ << ": "
                << matrix_.size() << " "
                << isRowSorted_ << " " << isColSorted_ 
                << std::endl;
*/
}


/*****Free Functions*MMreader*************************************************/
std::ostream& operator<<( std::ostream& os, std::tuple<int,int,double> data )
{
    os << "(" << std::get<0>(data) << "," << std::get<1>(data) << ") "
       << std::get<2>(data);

    return os;
}

template <typename T>
std::ostream& operator<<( std::ostream& os, std::vector<T>& vec )
{
    for (auto it=vec.begin(); it!=vec.end(); ++it)
    {
        os << *it << '\n';
    }

    os.flush();
    return os;
}

