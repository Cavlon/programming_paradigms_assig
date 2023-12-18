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

    for (size_t i = 1; i < static_cast<size_t>(argc); i++){
        int mod = (i-1) % dim;  // The component of the current vector

        // Remove the open square bracket of each new vector
        if (mod == 0){
            argv[i]++;
        } 

        float val = std::stof(argv[i]); // Cast each argument as a float

        m((i-1) / dim, mod) = val;  //(i-1)/dim is the column and mod is the row of the matrix
    }
    return m;
}