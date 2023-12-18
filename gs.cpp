#include "gs.h"

using namespace std;

// Get the Gram-Schmidt coefficient
float GetMu(const vector<float>& u, const vector<float>& v){
    float inner = u * v;

    // Avoid division by 0
    float norm = (inner != 0) ? v * v : 1;

    return inner / norm;
}

// Perform the Gram-Schmidt process without normalisation
void GS(const Matrix& m, Matrix& G, float coeffs[], size_t startInd){
    size_t dim = m.getDim();

    for (int i = startInd; i < dim; i++){
        G(i) = m(i);

        for (int j = 0; j < i; j++){
            float coeff = GetMu(m(i), G(j));
            coeffs[coeffInd(i, j)] = coeff;

            // Avoid unnecessary calculations if coefficient is zero
            if (coeff != 0.0f) {
                vector<float> proj = G(j) * coeff;

                G(i) = G(i) - proj;
            }
        }
    }
}