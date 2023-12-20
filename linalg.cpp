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
void Print(const vector<double>& v){
    for (double val : v){
        cout << val << " ";
    }
    cout << endl;
}

// Constructor using a specified dimension, all matrices are square
Matrix::Matrix(size_t d): dim(d), cols(d), vals(d, vector<double>(d, 0)){}

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
    vals.erase(vals.begin() + targetInd);

    vector<vector<double>> newVals;
    newVals.reserve(vals.size());

    for (size_t i = 0; i < vals.size(); ++i){
        newVals.push_back(vals[i]);
    }

    vals = move(newVals);

    --cols;
}

// Vector inner product
double operator*(const vector<double>& a, const vector<double>& b){
    // Numeric header has a built in dot product function, this initialises the result as 0
    return inner_product(a.begin(), a.end(), b.begin(), 0.0f);
}

// Vector subtraction
vector<double> operator-(const vector<double>& a, const vector<double>& b){
    vector<double> res(a.size());

    // Subtract from each element of a the corresponding element from b and store the result in res
    transform(a.begin(), a.end(), b.begin(), res.begin(), std::minus<double>());

    return res;
}

// Scalar vector multiplication
vector<double> operator*(const vector<double>& a, double b){
    vector<double> res(a.size());

    // Multiply each element of a by b and store it in res
    transform(a.begin(), a.end(), res.begin(), [b](double element) {
        return element * b;
    });

    return res;
}

void Insert(Matrix& m, size_t sourceInd, const size_t& targetInd){
    while (sourceInd > targetInd){
        swap(m(sourceInd), m(sourceInd-1));
        --sourceInd;
    }
}

bool IsNull(const vector<double>& v){
    for (const double& val : v){
        if (val != 0) return false;
    }
    return true;
}