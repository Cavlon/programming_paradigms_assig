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

        inline void Pop(size_t k){
            delete [] vals[k];
            vals.erase(vals.begin() + k);
            --cols;
        }

    private:
        size_t dim;
        size_t cols;
        std::vector<double*> vals;   // Vector of vectors storing the matrix values
};

inline double dot(const double* const & a, const double* const & b, const size_t& dim){
    double total = 0;
    for (size_t i = 0; i < dim; ++i){
        total += a[i] * b[i];
    }
    return total;
}

#endif