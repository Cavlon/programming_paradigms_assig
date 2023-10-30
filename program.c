#include <stdio.h>
#include <math.h>

void reduce(int* b1, int* b2, int dimension);
float norm(int* b, int dimension);
void swap_pointers(int** x, int** y);
void swap_nums(float *x, float *y);

int main(int argc, char** argv){
    int b1[] = {1, 1};
    int b2[] = {1, -1};

    b1[0] = *argv[1] - '0';
    b1[1] = *argv[2] - '0';
    b2[0] = *argv[3] - '0';
    b2[1] = *argv[4] - '0';

    reduce(&b1[0], &b2[0], 2);
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
    printf("%f\n", len1);
    printf("%f\n", len2);
    printf("%f\n", lennew);

    printf("\n%d ", b1[0]);
    printf("%d\n", b1[1]);
    printf("%d ", b2[0]);
    printf("%d\n\n", b2[1]);

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