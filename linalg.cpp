#include "linalg.h"

using namespace std;

// Constructor using a specified dimension, all matrices are square
Matrix::Matrix(int d): dim(d) {
    vals = new double*[d];
    for (int i = 0; i < d; ++i) {
        vals[i] = new double[d]();
    }
}

// Destructor frees memory from all vectors
Matrix::~Matrix() {
    for (int i = 0; i < dim; ++i) {
        delete [] vals[i];
    }
    delete [] vals;
}
