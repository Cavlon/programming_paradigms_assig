#include <cmath>
#include <string>
#include "linalg.h"

Matrix Parse(int argc, char** argv){
    int dim = (int)sqrt((double)(argc - 1));
    Matrix m(dim);

    int row = -1;
    for (int i=1; i<argc; i++){
        int mod = (i-1)%dim;

        if (mod == 0){
            argv[i]++;
            row++;
        } 
        float val = atof(argv[i]);
        m(row, mod) = val;
    }
    return m;
}