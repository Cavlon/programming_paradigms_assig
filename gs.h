#ifndef __GS
#define __GS
#include "linalg.h"

std::vector<float> Proj(const std::vector<float>& u, const std::vector<float>& v);
float GetMu(const std::vector<float>& u, const std::vector<float>& v);
Matrix GS(Matrix m);

#endif