#include <stdlib.h>
#include <stdio.h>
#include "gs.h"

mat create_matrix(int dim){
    mat X;
    X.vecs = malloc(dim * sizeof(vec));
    for (int i=0; i<dim; i++){
        X.vecs[i].vals = malloc(dim * sizeof(float));
    }
    return X;
}

float dot(float* u, float* v, int dim){
    float total = 0;
    for (int i=0; i<dim; i++){
        total += u[i] * v[i];
    }
    return total;
}

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

void copy_vec(float* u, float* v, int dim){
    for (int i=0; i<dim; i++){
        u[i] = v[i];
    }
}

void sub_vec(float* u, float* v, int dim){
    for (int i=0; i<dim; i++){
        u[i] -= v[i];
    }
}

float get_mu(mat X, mat G, int i, int j, int dim){
    vec u = X.vecs[i];
    vec v = G.vecs[j];
    return dot(u.vals, v.vals, dim) / dot(v.vals, v.vals, dim);
}

void swap_vecs(float* u, float* v, int dim){
    for (int i=0; i<dim; i++){
        float temp = u[i];
        u[i] = v[i];
        v[i] = temp;
    }
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

void print_vector(vec X, int dim){
    for (int i=0; i<dim; i++){
        printf("%f ", X.vals[i]);
    }
}

void print_matrix(mat X, int dim){
    for (int i=0; i<dim; i++){
        print_vector(X.vecs[i], dim);
        printf("\n");
    }
    printf("\n");
}

// int main(){
//     int dim = 3;
//     mat X = create_matrix(dim);

//     X.vecs[0].vals[0] = 1;
//     X.vecs[0].vals[1] = -1;
//     X.vecs[0].vals[2] = 1;
//     X.vecs[1].vals[0] = 1;
//     X.vecs[1].vals[1] = 0;
//     X.vecs[1].vals[2] = 1;
//     X.vecs[2].vals[0] = 1;
//     X.vecs[2].vals[1] = 1;
//     X.vecs[2].vals[2] = 2;

//     print_matrix(X, dim);

//     mat G = gs(X, dim);
//     print_matrix(G, dim);

//     return 0;
// }