#include <stdio.h>
#include <string.h>
#include <math.h>
#include "gs.h"

int dimension(int argc){
    return (int)sqrt((double)(argc - 1));
}

mat parse(int argc, char** argv, int dim){
    // printf("%d\n", argc);
    // printf("%d\n", dim);

    // for (int i=1; i<argc; i++){
    //     printf("%s ", argv[i]);
    // }

    // printf("\n");

    mat X = create_matrix(dim);
    int row = -1;
    for (int i=1; i<argc; i++){
        int mod = (i-1)%dim;

        if (mod == 0){
            argv[i]++;
            row++;
        } 
        float val = atof(argv[i]);
        X.vecs[row].vals[mod] = val;
        // printf("%f ", val);
    }
    // printf("\n\n");
    // print_matrix(X, dim);
    return X;
}

// int main(int argc, char** argv){
//     parse(argc, argv);
//     return 0;
// }