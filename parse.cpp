#include <cmath>
#include <string>
#include "linalg.h"

Matrix Parse(int argc, char** argv){
    int dim = static_cast<int>(sqrt(static_cast<double>(argc - 1)));

    // Check if the number of basis vectors = dimension of the basis i.e. a square matrix
    if (dim * dim != argc - 1) {
        throw std::runtime_error("The number of basis vectors is not equal to the dimension of the basis vectors");
    }

    // Create a square matrix of dimension m
    Matrix m(dim);

    size_t colInd = -1;

    double* col = nullptr;

    for (size_t i = 1; i < static_cast<size_t>(argc); ++i){
        size_t mod = (i-1) % dim;  // The component of the current vector

        if (mod == 0){
            if (*argv[i] != '['){
                throw std::runtime_error("Incorrect input format");
            }

            ++colInd;
            col = m(colInd);
            
            argv[i]++; // Remove the open square bracket of each new vector
        }

        double val = std::stod(argv[i]); // Cast each argument as a float

        col[mod] = val;  //(i-1)/dim is the column and mod is the row of the matrix
    }
    return m;
}