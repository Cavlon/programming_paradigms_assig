#ifndef __SVP
#define __SVP

#include <cmath>
#include "linalg.h"
#include "parse.h"

inline size_t coeffInd(int i, int j);
void LLL(Matrix& m, std::vector<double>& coeffs, std::vector<double>& norms, size_t& cols, const size_t& dim);
double Enumerate(const Matrix& m, const std::vector<double>& coeffs, const std::vector<double>& norms, const size_t& cols, const size_t& dim);
bool Init(const Matrix& m, std::vector<double>& coeffs, std::vector<double>& norms, const size_t& cols, const size_t& dim);
double SVP(Matrix& m);

#endif