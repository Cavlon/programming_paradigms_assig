#ifndef __LINALG
#define __LINALG

class Matrix
{
    public:
        Matrix(int d);   // Constructor
        ~Matrix();

        double*& operator()(int i){   // Vector retrieval
            return vals[i];
        }

        double* operator()(int i) const{   // Vector retrieval
            return vals[i];
        }

        inline int getDim() const{
            return dim;
        }  // Dimension getter

    private:
        int dim;
        double** vals;   // Vector of vectors storing the matrix values
};

inline double dot(const double* const & a, const double* const & b, const int& dim){
    double total = 0;
    for (int i = 0; i < dim; ++i){
        total += a[i] * b[i];
    }
    return total;
}

#endif