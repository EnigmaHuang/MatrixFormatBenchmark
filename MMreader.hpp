#ifndef MMREADER_HPP
#define MMREADER_HPP

#include <algorithm>
#include <tuple>
#include <vector>
#include <iostream>



template <typename T>
std::ostream& operator<<( std::ostream& os, std::vector<T>& vec );



/*****Class MMreader**********************************************************/
class MMreader
{
public:
    // Constructer reading a matrix market file
    MMreader(char const *fileName);

    // Geter
    bool isRowSorted() const { return isRowSorted_; }
    bool isColSorted() const { return isColSorted_; }
    void isRowSorted(bool status) { isRowSorted_ = status; }
    void isColSorted(bool status) { isColSorted_ = status; }
    int getRows() const { return M_; }
    int getCols() const { return N_; }
    int getNonZeros() const { return nz_; }
    std::vector< std::tuple<int,int,double> > const & getMatrx() const
    {
        return matrix_;
    }
    std::vector< std::tuple<int,int,double> >       & getMatrx()
    {
        return matrix_;
    }

private:
    int M_, N_, nz_;    // number of rows, collumns and nonzeros in Matrix
    bool isRowSorted_, isColSorted_;

    // intermidiate (cordinate based) representation of sparse matrix
    // zero based (!)
    std::vector< std::tuple<int,int,double> > matrix_;  

};


/*****Free Functions*MMreader*************************************************/
template <typename T>
std::ostream& operator<<( std::ostream& os, std::vector<T>& vec );

std::ostream& operator<<( std::ostream& os, std::tuple<int,int,double> data );

inline void sortByRow(MMreader& mmMatrix)
{
/*  std::cout << "sort by row" << std::endl;*/

    std::vector< std::tuple<int,int,double> > & matrix = mmMatrix.getMatrx(); 

    // first sort by cll
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

/*  std::cout << matrix;*/
}
#endif
