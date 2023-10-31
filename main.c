#include <stdio.h>
#include <string.h>
#include <math.h>
#include "gs.h"
#include "parse.h"

int main(int argc, char** argv){
    int dim = dimension(argc);

    mat X = parse(argc, argv, dim);
    print_matrix(X, dim);

    mat G = gs(X, dim);
    print_matrix(G, dim);
    return 0;
}