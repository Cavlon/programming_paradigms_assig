#include <cmath>
#include <string>
#include "linalg.h"

Matrix Parse(int argc, char** argv){
    int dim = static_cast<int>(sqrt(static_cast<double>(argc - 1)));

    // Check if the number of basis vectors = dimension of the basis i.e. a square matrix
    if (dim * dim != argc - 1) {
        std::cerr << "The number of basis vectors is not equal to the dimension of the basis vectors" << std::endl;
        std::exit(EXIT_FAILURE);
    }

    // Create a square matrix of dimension m
    Matrix m(dim);

    size_t colInd = -1;

    double* col = nullptr;

    for (size_t i = 1; i < static_cast<size_t>(argc); ++i){
        size_t mod = (i-1) % dim;  // The component of the current vector

        if (mod == 0){
            if (*argv[i] != '['){
                std::cerr << "Incorrect input format" << std::endl;
                std::exit(EXIT_FAILURE);
            }

            ++colInd;
            col = m(colInd);
            
            argv[i]++; // Remove the open square bracket of each new vector
        }


        double val = 0;
        try {
            val = std::stod(argv[i]); // Cast each argument as a float
        } catch (const std::invalid_argument& e) {
            std::cerr << "Input contains a non-number: " << argv[i] << std::endl;
            std::exit(EXIT_FAILURE);
        }

        col[mod] = val;
    }
    return m;
}