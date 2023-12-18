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
void Print(const vector<float>& v){
    for (float val : v){
        cout << val << " ";
    }
    cout << endl;
}

// Constructor using a specified dimension, all matrices are square
Matrix::Matrix(size_t d): dim(d), vals(d, vector<float>(d, 0)){}

// Getter for dimension member
size_t Matrix::getDim() const{
    return dim;
}

// Vector inner product
float operator*(const vector<float>& a, const vector<float>& b){
    // Numeric header has a built in dot product function, this initialises the result as 0
    return inner_product(a.begin(), a.end(), b.begin(), 0.0f);
}

// Vector subtraction
vector<float> operator-(const vector<float>& a, const vector<float>& b){
    vector<float> res(a.size());

    // Subtract from each element of a the corresponding element from b and store the result in res
    transform(a.begin(), a.end(), b.begin(), res.begin(), std::minus<float>());

    return res;
}

// Scalar vector multiplication
vector<float> operator*(const vector<float>& a, float b){
    vector<float> res(a.size());

    // Multiply each element of a by b and store it in res
    transform(a.begin(), a.end(), res.begin(), [b](float element) {
        return element * b;
    });

    return res;
}