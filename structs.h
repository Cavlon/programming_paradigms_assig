typedef struct vector
{
    float* vals;
} vec;

typedef struct matrix
{
    vec* vecs;
} mat;

mat create_matrix(int dim);
void print_vector(vec X, int dim);
void print_matrix(mat X, int dim);