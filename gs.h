#ifndef __GS
#define __GS
#include "linalg.h"

// Returns the index of the coefficient corresponding to m(i) and G(j)
inline size_t coeffInd(float i, int j) {
    return (((i - 1) * 0.5) * i) + j;
}
float GetMu(const std::vector<float>& u, const std::vector<float>& v);
void GS(const Matrix& m, Matrix& G, float coeffs[], size_t startInd = 0);

#endif