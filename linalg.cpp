#include "linalg.h"

using namespace std;

// Print a matrix
void Print(const Matrix& m) {
    int dim = m.getDim();
    for (size_t i = 0; i < dim; ++i) {
        for (size_t j = 0; j < dim; ++j) {
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

void Matrix::setCols(size_t newCols){
    cols = newCols;
}

void Matrix::Delete(size_t targetInd){

    size_t size = vals.size()-1;

    for (size_t i = targetInd; i < size; --i){
        vals.at(i) = vals.at(i+1);
    }

    delete [] vals[size];
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

double* add(const double* const & a, const double* const & b, const size_t& dim){
    double* res = new double[dim];
    transform(a, a+dim, b, res, std::plus<double>());
    return res;
}

double* scalar(const double* const & a, const double b, const size_t& dim){
    double* res = new double[dim];

    transform(a, a + dim, res, [b](double element) {
        return element * b;
    });

    return res;
}

void Insert(Matrix& m, size_t sourceInd, const size_t& targetInd){
    double* source = m(sourceInd);
    for (size_t i = sourceInd; i > targetInd; --i){
        m(i) = m(i - 1);
    }
    m(targetInd) = source;
}

bool IsNull(const double* const & v, const size_t& dim){
    for (size_t i = 0; i < dim; ++i){
        if (v[i] != 0) return false;
    }
    return true;
}