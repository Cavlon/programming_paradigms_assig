#include <iostream>
#include "linalg.h"

using namespace std;

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

void Print(const vector<float>& v){
    for (float val : v){
        cout << val << " ";
    }
    cout << endl;
}

Matrix::Matrix(size_t d): dim(d), vals(d, vector<float>(d, 0)){}

float operator*(const vector<float>& a, const vector<float>& b){
    float total = 0;
    for (int i = 0; i < a.size(); i++){
        total += a[i] * b[i];
    }
    return total;
}

float Matrix::getDim() const{
    return dim;
}

vector<float> operator-(const vector<float>& a, const vector<float>& b){
    vector<float> res(a.size());
    for (int i = 0; i < res.size(); i++){
        res[i] = a[i] - b[i];
    }
    return res;
}

vector<float> operator*(const vector<float>& a, float b){
    vector<float> res(a.size());
    for (int i = 0; i < res.size(); i++){
        res[i] = a[i] * b;
    }
    return res;
}

void Swap(vector<float>& a, vector<float>& b){
    vector<float> temp = a;
    a = b;
    b = temp;
}