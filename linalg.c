#include <stdio.h>
#include <stdlib.h>
#include "linalg.h"

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

void swap_vecs(float* u, float* v, int dim){
    for (int i=0; i<dim; i++){
        float temp = u[i];
        u[i] = v[i];
        v[i] = temp;
    }
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