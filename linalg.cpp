#include "linalg.h"

using namespace std;

// Constructor using a specified dimension, all matrices are square
Matrix::Matrix(size_t d): dim(d), cols(d){
    vals.reserve(d);
    for (size_t i = 0; i < d; ++i){
        vals.push_back(new double[d]);
    }
}

// Destructor frees memory from all vectors
Matrix::~Matrix() {
    for (size_t i = 0; i < cols; ++i){
        delete [] vals[i];
    }
}