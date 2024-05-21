#include <stdio.h>
#include "matrix.h"

float es[] = {
    1, 2, 
    3, 4,
}; 

int main(int argc, char const *argv[])
{
    Matrix A = __init__matrix(2, 2, es);
    Matrix B = __init__matrix(2, 2, es);
    // Matrix D = MAT_ID(2);
    // Matrix C = __alloc__matrix(2, 2);

    // __dot__matrix(C, A, B);
    __print__matrix(A);
    __print__matrix(B);
    __sum__matrix(A, B);
    __print__matrix(A);
    // __print__matrix(C);
    // __print__matrix(D);

    // Matrix E = MAT_MULTIPLY(3, A, B, C);
    // __print__matrix(E);
    return 0;
}
