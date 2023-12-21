#include <cmath>
#include <limits>
#include <chrono>
#include "linalg.h"
#include "parse.h"

using namespace std;

// Returns the index of the coefficient corresponding to m(i) and G(j)
inline size_t coeffInd(int i, int j) {
    return (static_cast<float>((i - 1) * 0.5) * i) + j;
}

void VectorSwap(Matrix& m, const vector<double>& norms, const vector<double>& coeffs, int& k, const size_t& dim){
    double c = dot(m(k), m(k), dim);
    size_t i = 0;

    while (i < k){
        // Print(m);
        // cout << "c = " << c << " norm = " << norms.at(i) << endl;
        // Print(coeffs);

        if (0.95 * norms.at(i) <= c){
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
    const size_t dim = m.getDim();
    unsigned int cols = m.getCols();
    vector<double> coeffs(coeffInd(cols, 0));  // Vector of Gram-Schmidt coefficients
    vector<double> norms(cols);  // Vector of the square norms of each Gram-Schmidt vector

    int k = 1;
    while (k < cols){

        // Initialisation of stage k
        norms.at(k) = dot(m(k), m(k), dim);
        if (k == 1) norms.at(0) = dot(m(0), m(0), dim);

        for (int j = 0; j < k; ++j){
            // cout << k << ' ' << j << endl;
            double s = dot(m(k), m(j), dim);
            
            for (int i = 0; i < j; ++i){
                s -= coeffs.at(coeffInd(j, i)) * coeffs.at(coeffInd(k, i)) * norms.at(i);
            }

            s = s / norms.at(j);

            norms.at(k) -= (s * s) * norms.at(j);
            coeffs.at(coeffInd(k, j)) = s;
        }
        // Print(norms);
        // Print(coeffs);
        // cout << endl;

        // Size reduction of the kth basis vector
        for (int j = k-1; j>=0; j--){
            double mu = coeffs.at(coeffInd(k, j));

            // Print(m(k), dim);
            // Print(m(j), dim);
            // cout << "k=" << k << " " << "j=" << j << " " << "mu=" << mu << '\n' << endl;
            

            if (abs(mu) > 0.5){
                mu = roundf(mu);            
                m(k) = sub(m(k), scalar(m(j), mu, dim), dim);
                for (int i = 0; i < j; ++i){
                    coeffs.at(coeffInd(k, i)) -= mu * coeffs.at(coeffInd(j, i));
                }
                coeffs.at(coeffInd(k, j)) -= mu;
            }
        }

        if (IsNull(m(k), dim)){
            m.Delete(k);
            k = 1;
            --cols;
            coeffs.resize(coeffInd(cols, 0), 0);
            norms.resize(cols, 0);
            continue;
        }

        // Basis vector swapping or incrementation of stage k
        VectorSwap(m, norms, coeffs, k, dim);
    }
}

chrono::microseconds SVP(int argc, char** argv){
    auto start = chrono::high_resolution_clock::now();

    Matrix m = Parse(argc, argv);
    // cout << "Input Matrix:\n";
    // Print(m);

    LLL(m);

    // cout << "Reduced Matrix:\n";
    // Print(m);

    float shortest = numeric_limits<float>::max();
    int shrtInd = 0;
    const size_t dim = m.getDim();
    const size_t cols = m.getCols();

    // Iterate through the reduced basis and find the smallest vector
    for (int i=0; i<cols; i++){
        float sqrnorm = dot(m(i), m(i), dim);
        if (sqrnorm < shortest){
            shortest = sqrnorm;
            shrtInd = i;
        }
    }

    // cout << "Shortest Vector is:\n";
    // Print(m(shrtInd), dim);
    // cout << "\nShortest Length: " << sqrtf(shortest) << endl;

    auto end = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::microseconds>(end - start);
    return duration;
}

int main(int argc, char** argv){

    size_t samples = 100000;
    chrono::microseconds time = SVP(argc, argv);

    for (size_t i = 0; i < samples-1; ++i){
        time += SVP(argc, argv);
    }

    cout << '\n' << "Time taken: " << time.count()/samples << " microseconds" << endl;
    
    return 0;
}