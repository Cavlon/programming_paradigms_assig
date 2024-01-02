#ifndef __SVP
#define __SVP

#include <iostream>
#include <fstream>
#include <cmath>
#include "linalg.h"
#include "parse.h"

inline size_t coeffInd(int i);
void LLL(Matrix* const m, std::vector<double>* const coeffs, std::vector<double>* const gNorms, std::vector<double>* const norms, size_t* const cols, const size_t& dim);
double Enumerate(Matrix* const m, std::vector<double>* const coeffs, std::vector<double>* const gNorms, const size_t& cols, const size_t& dim);
bool Init(Matrix* const m, std::vector<double>* const coeffs, std::vector<double>* const gNorms, std::vector<double>* const norms, const size_t& cols, const size_t& dim);
void SVP(Matrix* const m);

#endif