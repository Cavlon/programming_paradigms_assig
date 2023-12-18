#include <cmath>
#include <limits>
#include "linalg.h"
#include "gs.h"
#include "parse.h"

using namespace std;

// Performs LLL lattice reduction on the provided basis matrix
void LLL(Matrix& m){
    int dim = m.getDim();
    Matrix G(dim);  // The Gram-Schmidt matrix
    float coeffs[((dim - 1) >> 1) * dim ];  // Array of Gram-Schmidt coefficients

    GS(m, G, coeffs);

    int k = 1;
    while (k < dim){

        for (int j = k-1; j>=0; j--){
            float mu = coeffs[coeffInd(k, j)];

            // cout << "k=" << k << " " << "j=" << j << " " << "mu=" << mu << endl;
            // cout << coeffInd(k, j) << ' ' << coeffs[coeffInd(k, j)] << endl;
            // cout << GetMu(m(k), G(j)) << '\n' << endl;

            if (fabsf(mu) > 0.5){
                m(k) = m(k) - (m(j) * roundf(mu));
                GS(m, G, coeffs, k);
            }
        }

        vector<float> kcheck = G(k);
        vector<float> kchecklow = G(k-1);
        float mu = coeffs[coeffInd(k, k-1)];

        if ((kcheck * kcheck) > ((0.75 - (mu * mu)) * (kchecklow * kchecklow))){
            ++k;
        } else {
            swap(m(k), m(k-1));
            GS(m, G, coeffs, k-1);

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