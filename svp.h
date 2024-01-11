#ifndef __SVP
#define __SVP

#include <iostream>
#include <fstream>
#include <cmath>
#include "linalg.h"

inline size_t coeffInd(int i);
double Enumerate(Matrix* const m, double* const coeffs, double* const gNorms, const double* const vol, const size_t& dim);
void Init(Matrix* const m, double* const coeffs, double* const gNorms, double* const vol, const size_t& dim);
void SVP(Matrix* const m);

#endif