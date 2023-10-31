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

float dot(float* u, float* v, int dim);
vec proj(float* u, float* v, int dim);
void copy_vec(float* u, float* v, int dim);
void sub_vec(float* u, float* v, int dim);
mat gs(mat X, int dim);