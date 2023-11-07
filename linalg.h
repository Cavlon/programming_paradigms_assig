#ifndef __LINALG
#define __LINALG

typedef struct vector
{
    float* vals;
} vec;

typedef struct matrix
{
    vec* vecs;
} mat;

mat create_matrix(int dim);
float dot(float* u, float* v, int dim);
void copy_vec(float* u, float* v, int dim);
void sub_vec(float* u, float* v, int dim);
void swap_vecs(float* u, float* v, int dim);

void print_vector(vec X, int dim);
void print_matrix(mat X, int dim);

#endif