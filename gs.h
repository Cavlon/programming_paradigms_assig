#ifndef __GS
#define __GS
#include "linalg.h"

vec proj(float* u, float* v, int dim);
float get_mu(mat X, mat G, int i, int j, int dim);
mat gs(mat X, int dim);

#endif