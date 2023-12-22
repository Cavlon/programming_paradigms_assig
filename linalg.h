#ifndef __LINALG
#define __LINALG

#include <vector>
#include <algorithm> // Contains the transform function for operations on a range
#include <iostream>
#include <numeric> // Contains the inner product function

class Matrix
{
    public:
        Matrix(size_t d);   // Constructor
        ~Matrix();

        double& operator()(size_t i, size_t j)   // Value retrieval with 2 indices
        {
            return vals.at(i)[j];
        }

        double operator()(size_t i, size_t j) const  // Constant matrix value retrieval with 2 indices
        {
            return vals.at(i)[j];
        }

        double*& operator()(size_t i){   // Vector retrieval
            return vals.at(i);
        }

        double* operator()(size_t i) const{   // Vector retrieval
            return vals.at(i);
        }

        size_t getDim() const;  // Dimension getter

        size_t getCols() const;  // Dimension getter

        void setCols(size_t newCols);  // Dimension setter

        void Delete(size_t targetInd);

    private:
        size_t dim;
        size_t cols;
        std::vector<double*> vals;   // Vector of vectors storing the matrix values
};

void Print(const Matrix& m);
void Print(const double* const & v, size_t dim);
void Print(const std::vector<double>& v);
double dot(const double* const & a, const double* const & b, const size_t& dim);
double* add(const double* const & a, const double* const & b, const size_t& dim);
double* scalar(const double* const & a, const double b, const size_t& dim);
void Insert(Matrix& m, size_t sourceInd, const size_t& targetInd);
bool IsNull(const double* const & v, const size_t& dim);

#endif