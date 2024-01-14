#ifndef __SVP
#define __SVP

#include <iostream>
#include <fstream>
#include <cmath>
#include "linalg.h"

inline size_t coeffInd(int i);
double Enumerate(const Matrix& m, double* const coeffs, double* const gNorms, const double* const vol, const size_t& dim);
void Init(const Matrix& m, double* const coeffs, double* const gNorms, double* const vol, const size_t& dim);
void SVP(const Matrix& m);

#endif