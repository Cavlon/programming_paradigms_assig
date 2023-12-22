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

        if (0.99 * norms.at(i) <= c){
            double mu = coeffs.at(coeffInd(k, i));
            c -= (mu * mu) * norms.at(i);
            ++i;
        } else {
            Insert(m, k, i);
            k = --i;
            if (k < 1) k = 1;
            return;
        }
    }
    ++k;
}

// Performs LLL lattice reduction on the provided basis matrix
void LLL(Matrix& m, vector<double>& coeffs, vector<double>& norms, size_t& cols, const size_t& dim){

    int k = 1;
    while (k < cols){

        // Initialisation of stage k
        norms.at(k) = dot(m(k), m(k), dim);
        if (k == 1) norms.at(0) = dot(m(0), m(0), dim);

        for (int j = 0; j < k; ++j){
            double s = dot(m(k), m(j), dim);
            
            for (int i = 0; i < j; ++i){
                s -= coeffs.at(coeffInd(j, i)) * coeffs.at(coeffInd(k, i)) * norms.at(i);
            }

            s = s / norms.at(j);

            norms.at(k) -= (s * s) * norms.at(j);
            coeffs.at(coeffInd(k, j)) = s;
        }

        // Size reduction of the kth basis vector
        for (int j = k-1; j>=0; --j){
            double mu = coeffs.at(coeffInd(k, j));
            
            if (abs(mu) > 0.5){
                mu = roundf(mu);            
                m(k) = add(m(k), scalar(m(j), -mu, dim), dim);
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

void Enumerate(Matrix& m, const vector<double>& coeffs, const vector<double>& norms, const size_t& cols, const size_t& dim){

    double sqrRad = dot(m(0), m(0), dim);   // Radius for enumeration

    double* p = new double[cols+1]();   // Norms of the projections

    int* u = new int[cols]();   // Tentative integer coefficents for linear combination
    u[0] = 1;

    int* v = new int[cols]();   // Integer coefficents for linear combination
    v[0] = 1;

    double* c = new double[cols](); // Coefficients for linear combination

    unsigned int* w = new unsigned int[cols](); // Which direction to search through the level of the tree

    size_t k = 0;   // Level of the enumeration tree
    size_t nonZeroInd = 0;  // Last index that didn't have a 0 

    while (true){
        p[k] = p[k+1] + ((static_cast<double>(u[k]) - c[k]) * (static_cast<double>(u[k]) - c[k]) * norms.at(k));

        if (p[k] < sqrRad){
            if (k == 0) {
                sqrRad = p[k];
                for (size_t i = 0; i < cols; ++i){
                    v[i] = u[i];
                }

            } else {
                --k;

                c[k] = 0;
                for (size_t i = k + 1; i <= nonZeroInd; ++i){
                    c[k] -= coeffs.at(coeffInd(i, k)) * u[i];
                }

                u[k] = static_cast<int>(round(c[k]));
                w[k] = 1;
            }
        } else {
            ++k;

            if (k == cols) break;

            if (k >= nonZeroInd){
                nonZeroInd = k;
                u[k] += 1;
            } else {
                if (u[k] > c[k]) u[k] -= w[k];
                else u[k] += w[k];
                ++w[k];
            }
        }
    }
    double* res = new double[dim]();

    for (size_t i = 0; i < cols; ++i){
        res = add(res, scalar(m(i), v[i], dim), dim);
    }

    delete [] p;
    delete [] u;
    delete [] v;
    delete [] c;
    delete [] w;

    // cout << "Shortest Vector:\n";
    // Print(res, dim);

    delete [] res;
}

chrono::microseconds SVP(int argc, char** argv){
    auto start = chrono::high_resolution_clock::now();

    Matrix m = Parse(argc, argv);

    const size_t dim = m.getDim();
    size_t cols = m.getCols();
    vector<double> coeffs(coeffInd(cols, 0));  // Vector of Gram-Schmidt coefficients
    vector<double> norms(cols);  // Vector of the square norms of each Gram-Schmidt vector

    // LLL(m, coeffs, norms, cols, dim);

    for (size_t k = 1; k < cols; ++k){
        norms.at(k) = dot(m(k), m(k), dim);
        if (k == 1) norms.at(0) = dot(m(0), m(0), dim);

        for (int j = 0; j < k; ++j){
            double s = dot(m(k), m(j), dim);
            
            for (int i = 0; i < j; ++i){
                s -= coeffs.at(coeffInd(j, i)) * coeffs.at(coeffInd(k, i)) * norms.at(i);
            }

            s = s / norms.at(j);

            norms.at(k) -= (s * s) * norms.at(j);
            coeffs.at(coeffInd(k, j)) = s;
        }
    }

    // cout << "Reduced Matrix:\n";
    // Print(m);

    Enumerate(m, coeffs, norms, cols, dim);

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