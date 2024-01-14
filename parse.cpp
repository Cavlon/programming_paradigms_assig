#include "parse.h"

using namespace std;

Matrix Parse(int argc, char** argv) {
    // No basis given
    if (argc == 1) {
        cerr << "No input was given" << endl;
        exit(EXIT_FAILURE);
    }

    // 1 Dimensional Basis
    if (argc == 2) {
        if (*argv[1] != '[') {
            cerr << "Incorrect input format";
            exit(EXIT_FAILURE);
        }
        argv[1]++;

        double val = 0;
        try {
            val = stod(argv[1]);   // Cast each argument as a float
        } catch (const invalid_argument& e) {
            cerr << "Input contains a non-number: ";
            cerr << argv[1] << endl;
            exit(EXIT_FAILURE);
        }

        ofstream outFile("result.txt");

        if (!outFile.is_open()) {
            cerr << "Error opening the result file" << endl;
            outFile.close();
            exit(EXIT_FAILURE);
        }

        outFile << val;
        outFile.close();
        exit(EXIT_SUCCESS);
    }

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
                cerr << "Incorrect input format at: ";
                cerr << i-1 << endl;
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
