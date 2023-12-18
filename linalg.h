#ifndef __LINALG
#define __LINALG

#include <vector>

class Matrix
{
    public:
        Matrix(size_t d);

        float& operator()(size_t i, size_t j)
        {
            return vals[i][j];
        }

        float operator()(size_t i, size_t j) const
        {
            return vals[i][j];
        }

        std::vector<float>& operator()(size_t i){
            return vals[i];
        }

        float getDim() const;

        friend void Print(const Matrix& m);

    private:
        size_t dim;
        std::vector<std::vector<float>> vals;
};

void Print(const Matrix& m);
void Print(const std::vector<float>& v);
float operator*(const std::vector<float>& a, const std::vector<float>& b);
std::vector<float> operator-(const std::vector<float>& a, const std::vector<float>& b);
std::vector<float> operator*(const std::vector<float>& a, float b);
void Swap(std::vector<float>& a, std::vector<float>& b);

#endif