#include "svp.h"

using namespace std;

// Returns the index of the Gram-Schmidt coefficient of m(i) and m(0)
// All coefficients of m(i):  stored from coeffInd(i) to coeffInd(i) + (i - 1)
// Coefficent of m(i) & m(j) can be found with coeffInd(i) + j if j < i
inline size_t coeffInd(int i) {
    return (static_cast<float>((i - 1) * 0.5) * i);
}

double Enumerate(const Matrix& m, double* const coeffs,
                double* const gNorms, const double* const vol
                , const size_t& dim) {
    // Radius for enumeration
    double sqrRad = 0;

    // Norms of the projections
    double* p = new double[dim+1]();

    // Tentative integer coefficents for linear combination
    int* u = new int[dim]();
    u[0] = 1;

    // Integer coefficents for linear combination
    int* v = new int[dim]();
    v[0] = 0;

    // Coefficients for linear combination
    double* c = new double[dim]();

    // Direction to search through the tree
    unsigned int* w = new unsigned int[dim]();

    // Bound the search radius
    // Bound is an upper bound for the shortest norm
    // Pre-define low-dimensional constants
    if (dim < 11) {
        // Constants for 2-8 are Hermite coefficients
        // Constants for 9 & 10 are Hermite coefficient estimates
        const float* const vConstants = new float[9]{
            1.154701, 1.259922, 1.414214, 1.515717,
            1.665367, 1.811448, 2, 2.240647, 2.373267};

        sqrRad = pow(*vol, 1.0 / static_cast<double>(dim)) * vConstants[dim-2];

        delete [] vConstants;
    } else {
        // Calculate an estimate for the Hermite coefficient
        // calcConst = 2/PI
        constexpr double calcConst = 0.636619772368;

        double dimfp = static_cast<double>(dim);
        double gamRes = tgamma(2 + (dimfp / 2));
        gamRes = gamRes * gamRes * (*vol);

        sqrRad = calcConst * pow(gamRes, 1.0 / dimfp);
    }

    size_t k = 0;  // Level of the enumeration tree
    size_t nonZeroInd = 0;  // Last index that didn't have a 0

    while (true) {
        double diff = static_cast<double>(u[k]) - c[k];
        p[k] = p[k+1] +
        (diff * diff
        * gNorms[k]);

        if (p[k] < sqrRad) {
            if (k == 0) {   // New shortest vector found
                sqrRad = p[k];
                for (size_t i = 0; i < dim; ++i) {
                    v[i] = u[i];
                }

            } else {
                --k;

                c[k] = 0;
                // Minimise coeffInd calls
                // Index of first coeff related to vector k
                size_t ind = coeffInd(k+1) + k;
                for (size_t i = k + 1; i <= nonZeroInd; ++i) {
                    c[k] -= coeffs[ind] * u[i];
                    // Next index of coeff related to vector k
                    // It will always be stored i places after
                    // Where i is the current vector index
                    ind += i;
                }

                u[k] = static_cast<int>(round(c[k]));
                w[k] = 1;
            }
        } else {
            ++k;

            if (k == dim) break;

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
    for (size_t i = 0; i < dim*dim; ++i) {
        res[i % dim] += m.vals[i] * v[i/dim];
    }

    delete [] v;

    // Compute the shortest length from the shortest vector
    double total = sqrt(dot(res, res, dim));
    delete [] res;

    return total;
}

void Init(const Matrix& m, double* const coeffs,
            double* const gNorms, double* const vol,
            const size_t& dim) {
    for (int k = 0; k < dim; ++k) {
        double* kVec = &(m.vals[k*dim]);

        gNorms[k] = dot(kVec, kVec, dim);

        size_t kCoeffStart = coeffInd(k);
        for (int j = 0; j < k; ++j) {
            // Akin to the Gram-Schmidt process
            double s = dot(kVec, &(m.vals[j*dim]), dim);

            size_t jCoeffStart = coeffInd(j);
            for (int i = 0; i < j; ++i) {
                s -= coeffs[jCoeffStart + i] *
                coeffs[kCoeffStart + i] * gNorms[i];
            }

            s /= gNorms[j];

            // Calculate all  Gram-Schmidt norms
            gNorms[k] -= (s * s) * gNorms[j];

            // A linearly dependent vector exists
            if (gNorms[k] == 0) {
                cerr << "Linearly dependant, not a basis" << endl;
                exit(EXIT_FAILURE);
            }

            // Calculate all the Gram-Schmidt coefficients
            coeffs[kCoeffStart + j] = s;
        }

        // Calculate the squared volume
        *vol *= gNorms[k];
    }
}

void SVP(const Matrix& m) {
    const size_t dim = m.dim;

    // Gram-Schmidt coefficients
    double* coeffs = new double[coeffInd(dim)];

    // Square norms of each Gram-Schmidt vector
    double* gNorms = new double[dim];

    // The volume of the lattice's fundamental parallelepiped
    double* vol = new double;
    *vol = 1;

    Init(m, coeffs, gNorms, vol, dim);

    double res = Enumerate(m, coeffs, gNorms, vol, dim);

    delete vol;
    delete [] coeffs;
    delete [] gNorms;

    ofstream outFile("result.txt");

    if (!outFile.is_open()) {
        cerr << "Error opening the result file" << endl;
        outFile.close();
        exit(EXIT_FAILURE);
    }

    outFile << res;
    outFile.close();
}
