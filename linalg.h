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

        float& operator()(size_t i, size_t j)   // Value retrieval with 2 indices
        {
            return vals.at(i).at(j);
        }

        float operator()(size_t i, size_t j) const  // Constant matrix value retrieval with 2 indices
        {
            return vals.at(i).at(j);
        }

        std::vector<float>& operator()(size_t i){   // Vector retrieval
            return vals.at(i);
        }

        std::vector<float> operator()(size_t i) const{   // Constant matrix vector retrieval
            return vals.at(i);
        }

        size_t getDim() const;  // Dimension getter

        friend void Print(const Matrix& m); // Print matrix

    private:
        size_t dim;
        std::vector<std::vector<float>> vals;   // Vector of vectors storing the matrix values
};

void Print(const Matrix& m);
void Print(const std::vector<float>& v);
float operator*(const std::vector<float>& a, const std::vector<float>& b);
std::vector<float> operator-(const std::vector<float>& a, const std::vector<float>& b);
std::vector<float> operator*(const std::vector<float>& a, float b);

#endif