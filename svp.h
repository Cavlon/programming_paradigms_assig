#ifndef __SVP
#define __SVP

#include <cmath>
#include "linalg.h"
#include "parse.h"
#include <fstream>
#include <limits>

inline size_t coeffInd(int i);
void LLL(Matrix& m, std::vector<double>& coeffs, std::vector<double>& gNorms, std::vector<double>& norms, size_t& cols, const size_t& dim);
double Enumerate(const Matrix& m, const std::vector<double>& coeffs, const std::vector<double>& gNorms, const size_t& cols, const size_t& dim);
void SVP(Matrix& m);

#endif