#include "linalg.h"

using namespace std;

// Print a matrix
void Print(const Matrix& m) {
    int dim = m.getDim();
    int cols = m.getCols();
    for (size_t i = 0; i < dim; ++i) {
        for (size_t j = 0; j < cols; ++j) {
            cout << m(j,i) << " ";
        }
        cout << '\n';
    }
    cout << endl;
}

// Print a vector
void Print(const double* const & v, size_t dim){
    for (size_t i = 0; i < dim; ++i){
        cout << v[i] << " ";
    }
    cout << endl;
}

void Print(const vector<double>& v){
    for (double val : v){
        cout << val << " ";
    }
    cout << endl;
}

// Constructor using a specified dimension, all matrices are square
Matrix::Matrix(size_t d): dim(d), cols(d){
    vals.reserve(d);
    for (size_t i = 0; i < d; ++i){
        vals.push_back(new double[d]);
    }
}

Matrix::~Matrix() {
    for (size_t i = 0; i < cols; ++i){
        delete [] vals[i];
    }
}

// Getter for dimension member
size_t Matrix::getDim() const{
    return dim;
}

size_t Matrix::getCols() const{
    return cols;
}

void Matrix::Pop(){
    delete [] vals.at(cols - 1);
    vals.erase(vals.end() - 1);
    --cols;
}

double dot(const double* const & a, const double* const & b, const size_t& dim){
    double total = 0;
    for (size_t i = 0; i < dim; ++i){
        total += a[i] * b[i];
    }
    return total;
}

bool IsNull(const double* const & v, const size_t& dim){
    for (size_t i = 0; i < dim; ++i){
        if (v[i] != 0) return false;
    }
    return true;
}