#include <stdio.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include "linalg.h"
#include "gs.h"
#include "parse.h"

void LLL(mat X, int dim){
    mat G = gs(X, dim);
    // print_matrix(G, dim);

    int k = 1;
    while (k < dim){
        for (int j = k-1; j>=0; j--){
            float mu = get_mu(X, G, k, j, dim);
            // printf("k=%d, j=%d, mu=%f\n", k, j, mu);
            if (fabsf(mu) > 0.5){
                float* subvec = X.vecs[j].vals;
                for (int i=0; i<dim; i++){
                    float subval = subvec[i];
                    subval *= roundf(mu);
                    // printf("%f", subval);
                    // printf("\n");
                    X.vecs[k].vals[i] -= subval;
                }
                G = gs(X, dim);
            }
            // print_matrix(X, dim);
            // print_matrix(G, dim);
        }

        float* kcheck = G.vecs[k].vals;
        float* kchecklow = G.vecs[k-1].vals;
        float mu = get_mu(X, G, k, k-1, dim);
        if (dot(kcheck, kcheck, dim) > (0.75 - (mu * mu)) * dot(kchecklow, kchecklow, dim)){
            k++;
        } else {
            swap_vecs(X.vecs[k].vals, X.vecs[k-1].vals, dim);
            G = gs(X, dim);
            k = (k-1 > 1) ? k-1 : 1;
        }
    }
}

int main(int argc, char** argv){
    int dim = dimension(argc);

    mat X = parse(argc, argv, dim);
    printf("Transposed Input Matrix:\n");
    print_matrix(X, dim);

    printf("Transposed Reduced Matrix:\n");
    LLL(X, dim);
    print_matrix(X, dim);

    float shortest = FLT_MAX;
    int shrtInd = 0;
    for (int i=0; i<dim; i++){
        float* basisVec = X.vecs[i].vals;
        float sqrnorm = dot(basisVec, basisVec, dim);
        if (sqrnorm < shortest){
            shortest = sqrnorm;
            shrtInd = i;
        }
    }
    printf("Shortest Vector is:\n");
    print_vector(X.vecs[shrtInd], dim);
    printf("\nShortest Length: %f", sqrtf(shortest));
    return 0;
}