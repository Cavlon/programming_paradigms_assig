#include <stdlib.h>
#include "linalg.h"
#include "gs.h"

vec proj(float* u, float* v, int dim){
    float sqr_norm = dot(v, v, dim);
    float inner = dot(u, v, dim);

    float coeff = inner / sqr_norm;
    vec res;
    res.vals = malloc(dim * sizeof(float));

    for (int i=0; i<dim; i++){
        res.vals[i] = v[i] * coeff;
    }

    return res;
}

float get_mu(mat X, mat G, int i, int j, int dim){
    vec u = X.vecs[i];
    vec v = G.vecs[j];
    return dot(u.vals, v.vals, dim) / dot(v.vals, v.vals, dim);
}

mat gs(mat X, int dim){
    mat Y = create_matrix(dim);

    for (int i=0; i<dim; i++){
        copy_vec(Y.vecs[i].vals, X.vecs[i].vals, dim);
        for (int j=0; j<i; j++){
            vec proj_v = proj(X.vecs[i].vals, Y.vecs[j].vals, dim);
            sub_vec(Y.vecs[i].vals, proj_v.vals, dim);
        }
    }

    return Y;
}