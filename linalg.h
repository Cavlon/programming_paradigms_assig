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

        double& operator()(size_t i, size_t j)   // Value retrieval with 2 indices
        {
            return vals.at(i).at(j);
        }

        double operator()(size_t i, size_t j) const  // Constant matrix value retrieval with 2 indices
        {
            return vals.at(i).at(j);
        }

        std::vector<double>& operator()(size_t i){   // Vector retrieval
            return vals.at(i);
        }

        std::vector<double> operator()(size_t i) const{   // Constant matrix vector retrieval
            return vals.at(i);
        }

        size_t getDim() const;  // Dimension getter

        size_t getCols() const;  // Dimension getter

        void setCols(size_t newCols);  // Dimension setter

        void Delete(size_t targetInd);

    private:
        size_t dim;
        size_t cols;
        std::vector<std::vector<double>> vals;   // Vector of vectors storing the matrix values
};

void Print(const Matrix& m);
void Print(const std::vector<double>& v);
double operator*(const std::vector<double>& a, const std::vector<double>& b);
std::vector<double> operator-(const std::vector<double>& a, const std::vector<double>& b);
std::vector<double> operator*(const std::vector<double>& a, double b);
void Insert(Matrix& m, size_t sourceInd, const size_t& targetInd);
bool IsNull(const std::vector<double>& v);

#endif