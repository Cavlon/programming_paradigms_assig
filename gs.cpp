#include "gs.h"
#include <iostream>

using namespace std;

vector<float> Proj(const vector<float>& u, const vector<float>& v){
    float coeff = GetMu(u, v);
    return v * coeff;
}

float GetMu(const vector<float>& u, const vector<float>& v){
    float norm = v * v;
    float inner = u * v;
    return inner / norm;
}

Matrix GS(Matrix m){
    size_t dim = m.getDim();
    Matrix res(dim);

    for (int i = 0; i < dim; i++){
        res(i) = m(i);
        for (int j = 0; j < i; j++){
            vector<float> proj = Proj(m(i), res(j));
            res(i) = res(i) - proj;
        }
    }

    return res;
}