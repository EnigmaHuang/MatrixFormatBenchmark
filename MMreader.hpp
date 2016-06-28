#ifndef MMREADER_HPP
#define MMREADER_HPP

#include <tuple>
#include <vector>
#include <iostream>


/*****Class MMreader**********************************************************/
class MMreader
{
public:
    // Constructer reading a matrix market file
    MMreader(char const *fileName);

    // Geter
    bool isRowSorted()            const { return isRowSorted_; }
    bool isColSorted()            const { return isColSorted_; }
    void isRowSorted(bool status)       { isRowSorted_ = status; }
    void isColSorted(bool status)       { isColSorted_ = status; }
    int getRows()                 const { return M_; }
    int getCols()                 const { return N_; }
    int getNonZeros()             const { return matrix_.size(); }
    bool getSymmetry()            const { return isSymmetric_; }
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
    bool isSymmetric_;
    bool isRowSorted_, isColSorted_;

    // intermidiate (cordinate based) representation of sparse matrix
    // zero based (!)
    std::vector< std::tuple<int,int,double> > matrix_;

};




/*****Free Functions*MMreader*************************************************/
void sortByRow(MMreader& mmMatrix);

std::vector<int> getValsPerRow(MMreader const & mmMatrix);

std::vector< std::tuple<int,int> > getRowLengths(MMreader const & mmMatrix);

std::vector<int> getOffsets(std::vector<int> const & valuesPerLine);





/*****Output and Comperison Helper Funktions**********************************/
std::ostream& operator<<( std::ostream& os, std::tuple<int,int> data );
std::ostream& operator<<( std::ostream& os, std::tuple<int,int,double> data );

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

template <typename T>
bool operator ==(std::vector<T> const& a, T const& b)
{
    bool status = true;

    for (size_t i=0; i<a.size(); ++i)
    {
        status = (a[i] == b) && status;
    }

    return status;
}

template <typename T>
bool operator ==(T const& a, std::vector<T> const& b)
{
    return b==a;
}

template <typename T>
bool operator ==(std::vector<T> const& a, std::vector<T> const& b)
{
    bool status = a.size() == b.size();

    if (!status)
        return status;

    for (size_t i=0; i<a.size(); ++i)
    {
        status = (a[i] == b[i]) && status;
    }

    return status;
}

#endif
