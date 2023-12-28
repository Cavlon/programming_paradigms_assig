#ifndef __LINALG
#define __LINALG

#include <vector>
#include <iostream>
#include <cmath>

class Matrix
{
    public:
        Matrix(size_t d);   // Constructor
        ~Matrix();

        double& operator()(size_t i, size_t j)   // Value retrieval with 2 indices
        {
            return vals[i][j];
        }

        double operator()(size_t i, size_t j) const  // Constant matrix value retrieval with 2 indices
        {
            return vals[i][j];
        }

        double*& operator()(size_t i){   // Vector retrieval
            return vals[i];
        }

        double* operator()(size_t i) const{   // Vector retrieval
            return vals[i];
        }

        inline size_t getDim() const{
            return dim;
        }  // Dimension getter

        inline size_t getCols() const{
            return cols;
        };  // Column getter

        inline void Pop(){
            delete [] vals[cols - 1];
            vals.erase(vals.end() - 1);
            --cols;
        }

    private:
        size_t dim;
        size_t cols;
        std::vector<double*> vals;   // Vector of vectors storing the matrix values
};

void Print(const Matrix& m);
void Print(const double* const & v, size_t dim);
void Print(const std::vector<double>& v);

inline double dot(const double* const & a, const double* const & b, const size_t& dim){
    double total = 0;
    for (size_t i = 0; i < dim; ++i){
        total += a[i] * b[i];
    }
    return total;
}

inline bool IsNull(const double* const & v, const size_t& dim){
    for (size_t i = 0; i < dim; ++i){
        if (v[i] != 0) return false;
    }
    return true;
}

#endif