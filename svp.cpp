#include "svp.h"

using namespace std;

// Returns the index of the Gram-Schmidt coefficient of m(i) and m(0)
// All coefficients of m(i):  stored from coeffInd(i) to coeffInd(i) + (i - 1)
// Coefficent of m(i) & m(j) can be found with coeffInd(i) + j if j < i
inline size_t coeffInd(int i) {
    return (static_cast<float>((i - 1) * 0.5) * i);
}

// Performs LLL lattice reduction on the provided basis matrix
void LLL(Matrix* const m, vector<double>* const coeffs,
        vector<double>* const gNorms, vector<double>* const norms,
        size_t* const cols, const size_t& dim) {
    int k = 1;
    while (k < *cols) {
        // Initialisation of stage k
        double* kVec = (*m)(k);
        (*gNorms)[k] = (*norms)[k];
        if (k == 1) {
            (*gNorms)[0] = (*norms)[0];
        }

        size_t kCoeffStart = coeffInd(k);
        for (int j = 0; j < k; ++j) {
            double s = dot(kVec, (*m)(j), dim);

            size_t jCoeffStart = coeffInd(j);
            for (int i = 0; i < j; ++i) {
                s -= (*coeffs)[jCoeffStart + i] *
                (*coeffs)[kCoeffStart + i] * (*gNorms)[i];
            }

            s /= (*gNorms)[j];

            // Calculate all Gram-Schmidt norms
            (*gNorms)[k] -= (s * s) * (*gNorms)[j];
            // Calculate all the Gram-Schmidt coefficients
            (*coeffs)[kCoeffStart + j] = s;
        }

        // Size reduction of the kth basis vector
        bool changed = false;
        for (int j = k-1; j >= 0; --j) {
            double mu = (*coeffs)[kCoeffStart + j];

            // If true then the reduction calculations would be useless
            if (abs(mu) > 0.5) {
                changed = true;
                mu = round(mu);

                double* jVec = (*m)(j);
                for (size_t i = 0; i < dim; ++i) {
                    kVec[i] -= jVec[i] * mu;

                    // Rounds to prevent floating point errors
                    double roundVal = 1e-13;  // Rounding value
                    kVec[i] = round(kVec[i] / roundVal) * roundVal;
                }

                // Update all affected GS coefficients
                size_t jCoeffStart = coeffInd(j);
                for (int i = 0; i < j; ++i) {
                    (*coeffs)[kCoeffStart + i] -=
                    mu * (*coeffs)[jCoeffStart + i];
                }
                (*coeffs)[kCoeffStart + j] -= mu;
            }
        }
        if (changed) {  // Update the norm only if the vector has been changed
            (*norms)[k] = dot(kVec, kVec, dim);
        }

        if ((*norms)[k] == 0) {  // A zero vector has been found
            (*norms).erase((*norms).begin() + k);  // Remove the 0 vector's norm
            (*m).Pop(k);  // Remove the 0 vector

            // Resize the containers to account for the removal
            --(*cols);
            (*coeffs).resize(coeffInd((*cols)));
            (*gNorms).resize((*cols));
            continue;
        }

        // Lovász’ condition
        double checkCoeff = (*coeffs)[kCoeffStart + (k - 1)];
        // High delta ensures termination
        if ((0.999 - (checkCoeff * checkCoeff)) *
        (*gNorms)[k - 1] > (*gNorms)[k]) {
            swap((*m)(k), (*m)(k-1));
            swap((*norms)[k], (*norms)[k-1]);
            --k;
            if (k < 1) k = 1;
        } else {
            ++k;
        }
    }
}

double Enumerate(Matrix* const m, vector<double>* const coeffs,
                vector<double>* const gNorms,
                const size_t& cols, const size_t& dim) {
    // Radius for enumeration
    double sqrRad = (*gNorms)[0];

    // Norms of the projections
    double* p = new double[cols+1]();

    // Tentative integer coefficents for linear combination
    int* u = new int[cols]();
    u[0] = 1;

    // Integer coefficents for linear combination
    int* v = new int[cols]();
    v[0] = 1;

    // Coefficients for linear combination
    double* c = new double[cols]();

    // Direction to search through the tree
    unsigned int* w = new unsigned int[cols]();

    size_t k = 0;  // Level of the enumeration tree
    size_t nonZeroInd = 0;  // Last index that didn't have a 0

    while (true) {
        p[k] = p[k+1] +
        ((static_cast<double>(u[k]) - c[k]) * (static_cast<double>(u[k]) - c[k])
        * (*gNorms)[k]);

        if (p[k] < sqrRad) {
            if (k == 0) {   // New shortest vector found
                sqrRad = p[k];
                for (size_t i = 0; i < cols; ++i) {
                    v[i] = u[i];
                }

            } else {
                --k;

                c[k] = 0;
                for (size_t i = k + 1; i <= nonZeroInd; ++i) {
                    c[k] -= (*coeffs)[coeffInd(i) + k] * u[i];
                }

                u[k] = static_cast<int>(round(c[k]));
                w[k] = 1;
            }
        } else {
            ++k;

            if (k == cols) break;

            if (k >= nonZeroInd) {
                nonZeroInd = k;
                u[k] += 1;
            } else {
                if (u[k] > c[k]) {
                    u[k] -= w[k];
                } else {
                    u[k] += w[k];
                }
                ++w[k];
            }
        }
    }

    delete [] p;
    delete [] u;
    delete [] c;
    delete [] w;

    double* res = new double[dim]();

    // Compute the shortest vector as a linear combination
    for (size_t i = 0; i < cols; ++i) {
        double* iVec = (*m)(i);
        for (size_t j = 0; j < dim; ++j) {
            res[j] += iVec[j] * v[i];

            // Rounds to prevent floating point errors
            double roundVal = 1e-13;  // Rounding value
            res[j] = round(res[j] / roundVal) * roundVal;
        }
    }

    delete [] v;

    // Compute the shortest length from the shortest vector
    double total = 0;
    for (size_t i = 0; i < dim; ++i) {
        total += res[i] * res[i];
    }
    total = sqrt(total);

    return total;

    delete [] res;
}

bool Init(Matrix* const m, vector<double>* const coeffs,
            vector<double>* const gNorms, vector<double>* const norms,
            const size_t& cols, const size_t& dim) {
    for (int k = 0; k < cols; ++k) {
        double* kVec = (*m)(k);
        (*gNorms)[k] = (*norms)[k];
        if (k == 1) {
            (*gNorms)[0] = (*norms)[0];
        }

        size_t kCoeffStart = coeffInd(k);
        for (int j = 0; j < k; ++j) {
            double s = dot(kVec, (*m)(j), dim);

            size_t jCoeffStart = coeffInd(j);
            for (int i = 0; i < j; ++i) {
                s -= (*coeffs)[jCoeffStart + i] *
                (*coeffs)[kCoeffStart + i] * (*gNorms)[i];
            }

            s /= (*gNorms)[j];

            // Calculate all  Gram-Schmidt norms
            (*gNorms)[k] -= (s * s) * (*gNorms)[j];

            // A linearly dependent vector exists
            if ((*gNorms)[k] == 0) return true;

            // Calculate all the Gram-Schmidt coefficients
            (*coeffs)[kCoeffStart + j] = s;
        }
    }
    return false;
}

void SVP(Matrix* const m) {
    const size_t dim = (*m).getDim();
    size_t cols = (*m).getCols();
    // Gram-Schmidt coefficients
    vector<double> coeffs(coeffInd(cols));
    // Square norms of each Gram-Schmidt vector
    vector<double> gNorms(cols);
    // Square norms of each basis vector
    vector<double> norms(cols);

    // Pre-compute all the norms of the basis vectors
    for (size_t k = 0; k < cols; ++k) {
        double* kVec = (*m)(k);
        norms[k] = dot(kVec, kVec, dim);
    }

    if (Init(m, &coeffs, &gNorms, &norms, cols, dim)) {
        LLL(m, &coeffs, &gNorms, &norms, &cols, dim);   // Reduce the basis
    }

    double res = Enumerate(m, &coeffs, &gNorms, cols, dim);

    ofstream outFile("result.txt");

    if (!outFile.is_open()) {
        cerr << "Error opening the result file" << endl;
        outFile.close();
        exit(EXIT_FAILURE);
    }

    outFile << res;
    outFile.close();
}
