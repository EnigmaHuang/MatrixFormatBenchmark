
#include "MMreader.hpp"
#include "CSRMatrix.hpp"

#include <iostream>





int main(int argc, char *argv[])
{
    if (argc < 2)
	{
		fprintf(stderr, "Usage: %s [martix-market-filename]\n", argv[0]);
		exit(1);
	}

    MMreader test (argv[1]);
    CSR_Matrix csr_test(test);

    std::cout << "Matrix:"<< std::endl;
    std::cout << csr_test;


	return 0;
}

