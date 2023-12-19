#include <cmath>
#include <limits>
#include "linalg.h"
#include "parse.h"

using namespace std;

// Returns the index of the coefficient corresponding to m(i) and G(j)
inline size_t coeffInd(float i, int j) {
    return (((i - 1) * 0.5) * i) + j;
}

// Performs LLL lattice reduction on the provided basis matrix
void LLL(Matrix& m){
    unsigned int dim = m.getDim();
    float coeffs[((dim - 1) >> 1) * dim ];  // Array of Gram-Schmidt coefficients
    float norms[dim];

    int k = 1;
    while (k < dim){

        // Initialisation of stage k
        norms[k] = m(k) * m(k);
        if (k == 1) norms[0] = m(0) * m(0);

        for (int j = 0; j < k; ++j){
            float s = m(k) * m(j);
            
            for (int i = 0; i < j; ++i){
                s -= coeffs[coeffInd(j, i)] * coeffs[coeffInd(k, i)] * norms[i];
            }

            s = s / norms[j];

            norms[k] -= (s * s) * norms[j];

            coeffs[coeffInd(k, j)] = s;
        }

        // Size reduction of the kth basis vector
        for (int j = k-1; j>=0; j--){
            float mu = coeffs[coeffInd(k, j)];

            cout << "k=" << k << " " << "j=" << j << " " << "mu=" << mu << endl;
            // cout << coeffInd(k, j) << ' ' << coeffs[coeffInd(k, j)] << endl;

            if (fabsf(mu) > 0.5){
                mu = roundf(mu);
                m(k) = m(k) - (m(j) * mu);
                for (int i = 0; i < j; ++i){
                    coeffs[coeffInd(k, i)] -= mu * coeffs[coeffInd(j, i)];
                }
                coeffs[coeffInd(k, j)] -= mu;
            }
        }

        // Basis vector swapping or incrementation of stage k
        float mu = coeffs[coeffInd(k, k-1)];
        cout << "k norm = " << norms[k] << " k-1 norm = " << norms[k-1] << " mu = " << mu << " check = " << (0.75 - (mu * mu)) * norms[k-1] << '\n' << endl;

        if (norms[k] > ((0.75 - (mu * mu)) * norms[k-1])){
            ++k;
        } else {
            swap(m(k), m(k-1));

            --k;
            if (k < 1) k = 1;
        }
    }
}

int main(int argc, char** argv){
    Matrix m = Parse(argc, argv);
    cout << "Input Matrix:\n";
    Print(m);

    LLL(m);

    cout << "Reduced Matrix:\n";
    Print(m);

    float shortest = std::numeric_limits<float>::max();
    int shrtInd = 0;

    // Iterate through the reduced basis and find the smallest vector
    for (int i=0; i<m.getDim(); i++){
        const vector<float> basisVec = m(i);
        float sqrnorm = basisVec * basisVec;
        if (sqrnorm < shortest){
            shortest = sqrnorm;
            shrtInd = i;
        }
    }

    cout << "Shortest Vector is:\n";
    Print(m(shrtInd));
    cout << "\nShortest Length: " << sqrtf(shortest);
    return 0;
}