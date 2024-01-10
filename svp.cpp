#include "svp.h"

using namespace std;

// Returns the index of the Gram-Schmidt coefficient of m(i) and m(0)
// All coefficients of m(i):  stored from coeffInd(i) to coeffInd(i) + (i - 1)
// Coefficent of m(i) & m(j) can be found with coeffInd(i) + j if j < i
inline size_t coeffInd(int i) {
    return (static_cast<float>((i - 1) * 0.5) * i);
}

double Enumerate(Matrix* const m, double* const coeffs,
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

    // The constant to be multiplied with the volume
    // This constant is dimension based
    float constant = 1;

    // Pre-define low-dimensional constants
    if (dim < 11) {
        // Constants for 2-8 are Hermite coefficients
        // Constants for 9 & 10 are Hermite coefficient estimates
        const float* const vConstants = new float[9]{1.07457, 1.1224621, 1.189208, 1.231145, 1.290491, 1.34591, 1.414214, 1.496879, 1.540542};
        constant = vConstants[dim-2];   
        delete [] vConstants;
    } else {
        // Calculate an estimate for the Hermite coefficient
        // calcConst = root(2/PI)
        constexpr double calcConst = 0.797884560803;

        double dimfp = static_cast<double>(dim);
        constant = calcConst * pow(tgamma(2 + (dimfp / 2)), 1.0 / dimfp);
    }

    // Bound the search radius
    // Bound is an upper bound for the shortest norm
    sqrRad = constant * (*vol);
    sqrRad *= sqrRad;

    size_t k = 0;  // Level of the enumeration tree
    size_t nonZeroInd = 0;  // Last index that didn't have a 0

    while (true) {

        p[k] = p[k+1] +
        ((static_cast<double>(u[k]) - c[k]) * (static_cast<double>(u[k]) - c[k])
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
                for (size_t i = k + 1; i <= nonZeroInd; ++i) {
                    c[k] -= coeffs[coeffInd(i) + k] * u[i];
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
    for (size_t i = 0; i < dim; ++i) {
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
    double total = sqrt(dot(res, res, dim));
    delete [] res;

    return total;
}

void Init(Matrix* const m, double* const coeffs,
            double* const gNorms, double* const norms,
            double* const vol, const size_t& dim) {
    for (int k = 0; k < dim; ++k) {
        double* kVec = (*m)(k);

        size_t kCoeffStart = coeffInd(k);
        for (int j = 0; j < k; ++j) {
            double s = dot(kVec, (*m)(j), dim);

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

void SVP(Matrix* const m) {
    const size_t dim = (*m).getDim();
    // Gram-Schmidt coefficients
    double coeffs[coeffInd(dim)];
    // Square norms of each Gram-Schmidt vector
    double gNorms[dim];
    // Square norms of each basis vector
    double norms[dim];

    // Pre-compute all the norms of the basis vectors
    for (size_t k = 0; k < dim; ++k) {
        double* kVec = (*m)(k);
        double temp = dot(kVec, kVec, dim);
        norms[k] = temp;
        gNorms[k] = temp;
    }

    // The volume of the lattice's fundamental parallelepiped
    double* vol = new double;
    *vol = 1;

    Init(m, &coeffs[0], &gNorms[0], &norms[0], vol, dim);

    // Prepare volume for bound computation
    *vol = pow(*vol, 1.0 / static_cast<double>(dim*2));

    double res = Enumerate(m, &coeffs[0], &gNorms[0], vol, dim);

    delete vol;

    ofstream outFile("result.txt");

    if (!outFile.is_open()) {
        cerr << "Error opening the result file" << endl;
        outFile.close();
        exit(EXIT_FAILURE);
    }

    outFile << res;
    outFile.close();
}
