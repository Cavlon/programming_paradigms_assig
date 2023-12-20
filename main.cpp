#include <cmath>
#include <limits>
#include "linalg.h"
#include "parse.h"

using namespace std;

// Returns the index of the coefficient corresponding to m(i) and G(j)
inline size_t coeffInd(int i, int j) {
    return (static_cast<float>((i - 1) * 0.5) * i) + j;
}

void VectorSwap(Matrix& m, const vector<double>& norms, const vector<double>& coeffs, int& k){
    double c = m(k) * m(k);
    size_t i = 0;

    while (i < k){
        if (0.75 * norms.at(i) <= c){
            double mu = coeffs.at(coeffInd(k, i));
            c -= (mu * mu) * norms.at(i);
            ++i;
        } else {
            // cout << "Swap " << k << ' ' << i << '\n' << endl;
            Insert(m, k, i);
            k = --i;
            if (k < 1) k = 1;
            return;
        }
    }
    ++k;
}

// Performs LLL lattice reduction on the provided basis matrix
void LLL(Matrix& m){
    unsigned int cols = m.getCols();
    vector<double> coeffs(coeffInd(cols, 0));  // Vector of Gram-Schmidt coefficients
    vector<double> norms(cols);  // Vector of the square norms of each Gram-Schmidt vector

    int k = 1;
    while (k < cols){

        // Initialisation of stage k
        norms.at(k) = m(k) * m(k);
        if (k == 1) norms.at(0) = m(0) * m(0);

        for (int j = 0; j < k; ++j){
            // cout << k << ' ' << j << endl;
            double s = m(k) * m(j);
            
            for (int i = 0; i < j; ++i){
                s -= coeffs.at(coeffInd(j, i)) * coeffs.at(coeffInd(k, i)) * norms.at(i);
            }

            s = s / norms.at(j);

            norms.at(k) -= (s * s) * norms.at(j);
            coeffs.at(coeffInd(k, j)) = s;
        }

        // Size reduction of the kth basis vector
        for (int j = k-1; j>=0; j--){
            double mu = coeffs.at(coeffInd(k, j));

            // cout << "k=" << k << " " << "j=" << j << " " << "mu=" << mu << '\n' << endl;

            if (fabsf(mu) > 0.5){
                mu = roundf(mu);                
                m(k) = m(k) - (m(j) * mu);
                for (int i = 0; i < j; ++i){
                    coeffs.at(coeffInd(k, i)) -= mu * coeffs.at(coeffInd(j, i));
                }
                coeffs.at(coeffInd(k, j)) -= mu;
            }
        }

        if (IsNull(m(k))){
            m.Delete(k);
            k = 1;
            --cols;
            coeffs.resize(coeffInd(cols, 0), 0);
            norms.resize(cols, 0);
            continue;
        }

        // Basis vector swapping or incrementation of stage k
        VectorSwap(m, norms, coeffs, k);
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
    for (int i=0; i<m.getCols(); i++){
        const vector<double> basisVec = m(i);
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