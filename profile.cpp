#include <cmath>
#include "linalg.h"
#include "parse.h"

inline size_t coeffInd(int i);
void LLL(Matrix& m, std::vector<double>& coeffs, std::vector<double>& norms, size_t& cols, const size_t& dim);
double Enumerate(const Matrix& m, const std::vector<double>& coeffs, const std::vector<double>& norms, const size_t& cols, const size_t& dim);
bool Init(const Matrix& m, std::vector<double>& coeffs, std::vector<double>& norms, const size_t& cols, const size_t& dim);
double SVP(Matrix& m);

using namespace std;

// Returns the index of the coefficient corresponding to m(i) and G(j)
inline size_t coeffInd(int i) {
    return (static_cast<float>((i - 1) * 0.5) * i);
}

// Performs LLL lattice reduction on the provided basis matrix
void LLL(Matrix& m, vector<double>& coeffs, vector<double>& norms, size_t& cols, const size_t& dim){

    int k = 1;
    while (k < cols){

        // Initialisation of stage k
        double* kVec = m(k);
        norms[k] = dot(kVec, kVec, dim);
        if (k == 1){
            double* firstVec = m(0);
            norms[0] = dot(firstVec, firstVec, dim);
        }

        size_t kCoeffStart = coeffInd(k);
        for (int j = 0; j < k; ++j){
            double s = dot(kVec, m(j), dim);
            
            size_t jCoeffStart = coeffInd(j);
            for (int i = 0; i < j; ++i){
                s -= coeffs[jCoeffStart + i] * coeffs[kCoeffStart + i] * norms[i];
            }

            s /= norms[j];

            norms[k] -= (s * s) * norms[j];
            coeffs[kCoeffStart + j] = s;
        }

        // Size reduction of the kth basis vector
        for (int j = k-1; j>=0; --j){
            double mu = coeffs[kCoeffStart + j];
            
            if (abs(mu) > 0.5){
                mu = round(mu);

                double* jVec = m(j);
                for (size_t i = 0; i < dim; ++i){
                    kVec[i] -= jVec[i] * mu;

                    double roundVal = 1e-13;  // Rounding value
                    kVec[i] = round(kVec[i] / roundVal) * roundVal;
                }

                size_t jCoeffStart = coeffInd(j);
                for (int i = 0; i < j; ++i){
                    coeffs[kCoeffStart + i] -= mu * coeffs[jCoeffStart + i];
                }
                coeffs[kCoeffStart + j] -= mu;
            }
        }

        if (IsNull(kVec, dim)){
            swap(m(k), m(cols-1));
            m.Pop();
            k = 1;
            --cols;
            coeffs.resize(coeffInd(cols), 0);
            norms.resize(cols, 0);
            continue;
        }

        // Basis vector swapping or incrementation of stage k
        size_t checkCoeff = coeffs[kCoeffStart + (k - 1)];
        if ((0.99 - (checkCoeff * checkCoeff)) * norms[k - 1] > norms[k]){
            swap(m(k), m(k-1));
            --k;
            if (k < 1) k = 1;
        } else ++k;
    }
}

double Enumerate(const Matrix& m, const vector<double>& coeffs, const vector<double>& norms, const size_t& cols, const size_t& dim){

    double sqrRad = norms[0];   // Radius for enumeration

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
        p[k] = p[k+1] + ((static_cast<double>(u[k]) - c[k]) * (static_cast<double>(u[k]) - c[k]) * norms[k]);

        if (p[k] < sqrRad ){
            if (k == 0) {
                sqrRad = p[k];
                for (size_t i = 0; i < cols; ++i){
                    v[i] = u[i];
                }

            } else {
                --k;

                c[k] = 0;
                for (size_t i = k + 1; i <= nonZeroInd; ++i){
                    c[k] -= coeffs[coeffInd(i) + k] * u[i];
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

        double* iVec = m(i);
        for (size_t j = 0; j < dim; ++j){
            res[j] += iVec[j] * v[i];

            double roundVal = 1e-13;  // Rounding value
            res[j] = round(res[j] / roundVal) * roundVal;
        }
    }

    delete [] p;
    delete [] u;
    delete [] v;
    delete [] c;
    delete [] w;

    double total = 0;
    for (size_t i = 0; i < dim; ++i){
        total += res[i] * res[i];
    }
    total = sqrt(total);

    return total;

    delete [] res;
}

bool Init(const Matrix& m, vector<double>& coeffs, vector<double>& norms, const size_t& cols, const size_t& dim){
    size_t kCoeffStart = 0;
    size_t jCoeffStart = 0;

    for (size_t k = 1; k < cols; ++k){
        double* kVec = m(k);
        norms[k] = dot(kVec, kVec, dim);
        if (k == 1){
            double* firstVec = m(0);
            norms[0] = dot(firstVec, firstVec, dim);
        }

        kCoeffStart = coeffInd(k);
        for (int j = 0; j < k; ++j){
            double s = dot(kVec, m(j), dim);
            
            jCoeffStart = coeffInd(j);
            for (int i = 0; i < j; ++i){
                s -= coeffs[jCoeffStart + i] * coeffs[kCoeffStart + i] * norms[i];
            }

            s = s / norms[j];

            norms[k] -= (s * s) * norms[j];
            if (norms[k] == 0) return true;  // A vector is linearly dependent

            coeffs[kCoeffStart + j] = s;
        }
    }
    return false;
}

double SVP(Matrix& m){

    const size_t dim = m.getDim();
    size_t cols = m.getCols();
    vector<double> coeffs(coeffInd(cols));  // Vector of Gram-Schmidt coefficients
    vector<double> norms(cols);  // Vector of the square norms of each Gram-Schmidt vector

    bool ld = Init(m, coeffs, norms, cols, dim);

    if (ld){
        LLL(m, coeffs, norms, cols, dim);
    }

    return Enumerate(m, coeffs, norms, cols, dim);
}

int main(int argc, char** argv){

    size_t samples = 100000;

    for (size_t i = 0; i < samples; ++i){
        Matrix m = Parse(argc, argv);
        SVP(m);
    }
    
    return 0;
}