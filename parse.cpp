#include "parse.h"

using namespace std;

Matrix Parse(int argc, char** argv) {
    int dim = static_cast<int>(sqrt(static_cast<double>(argc - 1)));

    // Check if the number of basis vectors = dimension
    if (dim * dim != argc - 1) {
        cerr << "Input is not a square basis" << endl;
        exit(EXIT_FAILURE);
    }

    // Create a square matrix of dimension m
    Matrix m(dim);

    size_t colInd = -1;

    double* col = nullptr;

    for (size_t i = 1; i < static_cast<size_t>(argc); ++i) {
        size_t mod = (i-1) % dim;  // The component of the current vector

        if (mod == 0) {
            if (*argv[i] != '[') {
                cerr << "Incorrect input format" << endl;
                exit(EXIT_FAILURE);
            }

            ++colInd;
            col = m(colInd);

            argv[i]++;  // Remove the open square bracket of each new vector
        }

        double val = 0;
        try {
            val = stod(argv[i]);   // Cast each argument as a float
        } catch (const invalid_argument& e) {
            cerr << "Input contains a non-number: ";
            cerr << argv[i] << endl;
            exit(EXIT_FAILURE);
        }

        col[mod] = val;
    }

    return m;
}
