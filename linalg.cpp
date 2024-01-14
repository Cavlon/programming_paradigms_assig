#include "linalg.h"

using namespace std;

// Constructor using a specified dimension, all matrices are square
Matrix::Matrix(int d): dim(d), vals(new double[d*d]()) {}

// Destructor frees memory from the matrix
Matrix::~Matrix() {
    delete [] vals;
}
