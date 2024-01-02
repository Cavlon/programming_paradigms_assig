#include "svp.h"

using namespace std;

// Returns the index of the Gram-Schmidt coefficient of m(i) and m(0)
// All coefficients of a vector m(i) are stored from coeffInd(i) to coeffInd(i) + (i - 1) in the coeffs vector
// The coefficent between m(i) and m(j) can be found with coeffInd(i) + j if j < i
inline size_t coeffInd(int i) {
    return (static_cast<float>((i - 1) * 0.5) * i);
}

// Performs LLL lattice reduction on the provided basis matrix
void LLL(Matrix& m, vector<double>& coeffs, vector<double>& gNorms, vector<double>& norms, size_t& cols, const size_t& dim){

    int k = 1;
    while (k < cols){

        // Initialisation of stage k
        double* kVec = m(k);
        gNorms[k] = norms[k];
        if (k == 1){
            gNorms[0] = norms[0];
        }

        size_t kCoeffStart = coeffInd(k);
        for (int j = 0; j < k; ++j){
            double s = dot(kVec, m(j), dim);
            
            size_t jCoeffStart = coeffInd(j);
            for (int i = 0; i < j; ++i){
                s -= coeffs[jCoeffStart + i] * coeffs[kCoeffStart + i] * gNorms[i];
            }

            s /= gNorms[j];

            gNorms[k] -= (s * s) * gNorms[j];   // Calculate the Gram-Schmidt norm of the kth vector
            coeffs[kCoeffStart + j] = s;    // Calculate all the Gram-Schmidt coefficients of the kth vector with each preceeding vector
        }

        // Size reduction of the kth basis vector
        bool changed = false;
        for (int j = k-1; j>=0; --j){
            double mu = coeffs[kCoeffStart + j];
            
            if (abs(mu) > 0.5){ // If the coefficient would round to 0 then the reduction calculations would be useless
                changed = true;
                mu = round(mu);

                double* jVec = m(j);
                for (size_t i = 0; i < dim; ++i){
                    kVec[i] -= jVec[i] * mu;

                    double roundVal = 1e-13;  // Rounding value
                    kVec[i] = round(kVec[i] / roundVal) * roundVal; // Rounds the value slightly to prevent floating point errors
                }

                size_t jCoeffStart = coeffInd(j);
                for (int i = 0; i < j; ++i){
                    coeffs[kCoeffStart + i] -= mu * coeffs[jCoeffStart + i];    // Update all affected GS coefficients
                }
                coeffs[kCoeffStart + j] -= mu;  // Update all affected GS coefficients
            }
        }
        if (changed){   // Update the norm only if the vector has been changed
            norms[k] = dot(kVec, kVec, dim);
        }

        if (norms[k] == 0){ // A zero vector has been found (either a 0 vector in the basis or a linearly dependant vector was found)
            norms.erase(norms.begin() + k); // Remove the 0 vector's norm
            m.Pop(k);   // Remove the 0 vector

            // Resize the containers to account for the removal
            --cols;
            coeffs.resize(coeffInd(cols));
            gNorms.resize(cols);
            continue;
        }

        // Basis vector swapping or incrementation of stage k (Lovász’ condition)
        double checkCoeff = coeffs[kCoeffStart + (k - 1)];
        if ((0.999 - (checkCoeff * checkCoeff)) * gNorms[k - 1] > gNorms[k]){    // High delta chosen to ensure successful termination
            swap(m(k), m(k-1));
            swap(norms[k], norms[k-1]);
            --k;
            if (k < 1) k = 1;
        } else ++k;
    }
}

double Enumerate(const Matrix& m, const vector<double>& coeffs, const vector<double>& gNorms, const size_t& cols, const size_t& dim){

    double sqrRad = gNorms[0];   // Radius for enumeration

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
        p[k] = p[k+1] + ((static_cast<double>(u[k]) - c[k]) * (static_cast<double>(u[k]) - c[k]) * gNorms[k]);

        if (p[k] < sqrRad ){
            if (k == 0) {   // New shortest vector found
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

    delete [] p;
    delete [] u;
    delete [] c;
    delete [] w;

    double* res = new double[dim]();

    for (size_t i = 0; i < cols; ++i){  // Compute the shortest vector as a linear combination

        double* iVec = m(i);
        for (size_t j = 0; j < dim; ++j){
            res[j] += iVec[j] * v[i];

            double roundVal = 1e-13;  // Rounding value
            res[j] = round(res[j] / roundVal) * roundVal;   // Rounds the value slightly to prevent floating point errors
        }
    }

    delete [] v;

    double total = 0;
    for (size_t i = 0; i < dim; ++i){   // Compute the shortest length from the shortest vector
        total += res[i] * res[i];
    }
    total = sqrt(total);

    return total;

    delete [] res;
}

bool Init(const Matrix& m, vector<double>& coeffs, vector<double>& gNorms, vector<double>& norms, const size_t& cols, const size_t& dim){
    for (int k = 0; k < cols; ++k){
        double* kVec = m(k);
        gNorms[k] = norms[k];
        if (k == 1){
            gNorms[0] = norms[0];
        }

        size_t kCoeffStart = coeffInd(k);
        for (int j = 0; j < k; ++j){
            double s = dot(kVec, m(j), dim);
            
            size_t jCoeffStart = coeffInd(j);
            for (int i = 0; i < j; ++i){
                s -= coeffs[jCoeffStart + i] * coeffs[kCoeffStart + i] * gNorms[i];
            }

            s /= gNorms[j];

            gNorms[k] -= (s * s) * gNorms[j];   // Calculate the Gram-Schmidt norm of the kth vector

            if (gNorms[k] == 0) return true;
            
            coeffs[kCoeffStart + j] = s;    // Calculate all the Gram-Schmidt coefficients of the kth vector with each preceeding vector
        }
    }
    return false;
}

void SVP(Matrix& m){

    const size_t dim = m.getDim();
    size_t cols = m.getCols();
    vector<double> coeffs(coeffInd(cols));  // Vector of Gram-Schmidt coefficients
    vector<double> gNorms(cols);  // Vector of the square norms of each Gram-Schmidt vector
    vector<double> norms(cols); // Vector of the square norms of each basis vector

    for (size_t k = 0; k < cols; ++k){  // Pre-compute all the norms of the basis vectors
        double* kVec = m(k);
        norms[k] = dot(kVec, kVec, dim);
    }

    if (Init(m, coeffs, gNorms, norms, cols, dim)){
        LLL(m, coeffs, gNorms, norms, cols, dim);   // Reduce the basis
    }

    double res = Enumerate(m, coeffs, gNorms, cols, dim);

    ofstream outFile("result.txt");

    if (!outFile.is_open()){
        cerr << "Error opening the result file" << endl;
        outFile.close();
        exit(EXIT_FAILURE);
    }

    outFile << res;
    outFile.close();
}