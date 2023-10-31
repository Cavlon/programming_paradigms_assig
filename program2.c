#include <stdio.h>
#include <math.h>
#include "structs.h"

void reduce(int* b1, int* b2, int dimension);
float norm(int* b, int dimension);
void swap_pointers(int** x, int** y);
void swap_nums(float *x, float *y);
void print_vector(int* vec, int dimension);

int main(int argc, char** argv){
    int b1[] = {2, 1, 1};
    int b2[] = {1, -4, 1};
    int b3[] = {-6, -1, 2};

    // b1[0] = *argv[1] - '0';
    // b1[1] = *argv[2] - '0';
    // b1[2] = *argv[3] - '0';
    // b2[0] = *argv[4] - '0';
    // b2[1] = *argv[5] - '0';
    // b2[2] = *argv[6] - '0';
    // b3[0] = *argv[7] - '0';
    // b3[1] = *argv[8] - '0';
    // b3[2] = *argv[9] - '0';

    reduce(&b1[0], &b2[0], 3);
    reduce(&b1[0], &b3[0], 3);
}

void reduce(int* b1, int* b2, int dimension){
    int new[dimension];

    float len1 = norm(b1, dimension);
    float len2 = norm(b2, dimension);

    if (len2 > len1) {
        swap_pointers(&b1, &b2);
        swap_nums(&len1, &len2);
    }

    for (int i=0; i<dimension; i++){
        new[i] = b1[i] - b2[i];
    }

    float lennew = norm(&new[0], dimension);
    printf("\n%f\n", len1);
    printf("%f\n", len2);
    printf("%f\n", lennew);

    print_vector(b1, dimension);
    print_vector(b2, dimension);

    if (lennew > len1) return;

    for (int i=0; i<dimension; i++){
        b1[i] = new[i];
    }

    reduce(b2, b1, dimension);
}

float norm(int* b, int dimension){
    float total = 0;
    for (int i=0; i<dimension; i++){
        total += b[i] * b[i];
    }
    return sqrtf(total);
}

void swap_pointers(int **x, int **y) {
  int *tmp;
  tmp = *x;
  *x = *y;
  *y = tmp;
}

void swap_nums(float *x, float *y) {
  float tmp;
  tmp = *x;
  *x = *y;
  *y = tmp;
}

void print_vector(int* vec, int dimension){
    printf("\n");
    for (int i=0; i<dimension; i++){
        printf("%d ", vec[i]);
    }
}