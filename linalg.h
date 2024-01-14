#ifndef __LINALG
#define __LINALG

class Matrix
{
    public:
        Matrix(int d);   // Constructor
        ~Matrix();

        // Dimension of the matrix
        int dim;
        // 1D array of the matrix values
        // Stored sequentially as the columns
        // i.e. [vec1 ..., vec2 ..., vec3 ..., etc.]
        // Index of 1st component of ith vector = i*dim
        double* vals;
};

// The dot product of 2 vectors
inline double dot(const double* const & a, const double* const & b, const int& dim){
    double total = 0;
    for (int i = 0; i < dim; ++i){
        total += a[i] * b[i];
    }
    return total;
}

#endif