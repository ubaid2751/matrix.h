#ifndef matrix_h
#define matrix_h

#include <stdlib.h>
#include <stddef.h>
#include <math.h>
#include <stdio.h>
#include <stdarg.h>
#include "activations.h"

#ifndef MAT_ASSERT
#include <assert.h>
#define MAT_ASSERT assert
#endif

#define __NEWLINE__ printf("\n") 

//Defined Matrix
typedef struct {
    size_t rows;
    size_t cols;
    size_t stride;
    float* elements;
} Matrix;

// Defined Macros for Matrix
#define MAT_AT(mat, row, col) (mat).elements[(row)*(mat).stride + (col)]
#define __print__matrix(mat) __display__matrix(mat, #mat)
#define MAT_ID(size) __identity__matrix(size)
#define MAT_CONV(AB) __convolution__matrix(A, B)
#define MAT_MULTIPLY(num_args, ...) __multiply__matrix(num_args, __VA_ARGS__)

// Functions with return type as Matrix
Matrix __alloc__matrix(size_t rows, size_t cols);
Matrix __init__matrix(size_t rows, size_t cols, float* elements);
Matrix __init__random__matrix(size_t rows, size_t cols);
Matrix __identity__matrix(size_t size);
Matrix __convolution__matrix(Matrix A, Matrix B);
Matrix __multiply__matrix(int num_args, ...);

// Functions with no return type
void __display__matrix(Matrix mat, const char* name);
void __sum__matrix(Matrix dest, Matrix mat);
static void __dot__matrix(Matrix dest, Matrix a1, Matrix a2); 
void __copy__matrix(Matrix dest, Matrix src);
void __sigmoidf__matrix(Matrix mat);
void __relu__matrix(Matrix mat);
void __scalar__multiply(Matrix mat, float x);

// Functions with return type as float
float _random_float_();

// Functions written:
float _random_float_() {
    return (float) rand() / (float) RAND_MAX;
}

Matrix __alloc__matrix(size_t rows, size_t cols) {
    Matrix mat;
    mat.rows = rows;
    mat.cols = cols;
    mat.stride = cols;
    mat.elements = (float*)malloc(rows * cols * sizeof(float));

    MAT_ASSERT(mat.elements != NULL);

    for (size_t row = 0; row < mat.rows; row++) {
        for (size_t col = 0; col < mat.cols; col++) {
            MAT_AT(mat, row, col) = 0.0;
        }
    }

    return mat;
}

Matrix __init__matrix(size_t rows, size_t cols, float* elements) {
    Matrix mat = __alloc__matrix(rows, cols);
    if (elements) {
        for (size_t i = 0; i < rows * cols; i++) {
            mat.elements[i] = elements[i];
        }
    }
    return mat;
}

Matrix __init__random__matrix(size_t rows, size_t cols) {
    Matrix mat = __alloc__matrix(rows, cols);
    for (size_t row = 0; row < mat.rows; row++) {
        for (size_t col = 0; col < mat.cols; col++) {
            MAT_AT(mat, row, col) = _random_float_();
        }
    }
    return mat;
}

Matrix __identity__matrix(size_t size) {
    Matrix mat = __alloc__matrix(size, size);
    for (size_t i = 0; i < size; i++) {
        MAT_AT(mat, i, i) = 1.0;
    }
    return mat;
}

Matrix __convolution__matrix(Matrix A, Matrix B) {
    size_t rowC = A.rows - B.rows + 1;
    size_t colC = A.cols - B.cols + 1;
    Matrix C = __alloc__matrix(rowC, colC);

    for (size_t row = 0; row < rowC; row++) {
        for (size_t col = 0; col < colC; col++) {
            float sum = 0.0;
            for (size_t i = 0; i < B.rows; i++) {
                for (size_t j = 0; j < B.cols; j++) {
                    sum += MAT_AT(A, row + i, col + j) * MAT_AT(B, i, j);
                }
            }
            MAT_AT(C, row, col) = sum;    
        }
    }
    return C;
}

void __display__matrix(Matrix mat, const char* name) {
    printf("%s = {\n", name);
    for (size_t row = 0; row < mat.rows; row++) {
        for (size_t col = 0; col < mat.cols; col++) {
            printf("\t%f ", MAT_AT(mat, row, col)); 
        }
        __NEWLINE__;
    }
    printf("}");
    __NEWLINE__;
}

void __sum__matrix(Matrix dest, Matrix mat) {
    for (size_t row = 0; row < mat.rows; row++) {
        for (size_t col = 0; col < mat.cols; col++) {
            MAT_AT(dest, row, col) += MAT_AT(mat, row, col);
        }
    }
}

static void __dot__matrix(Matrix dest, Matrix a1, Matrix a2) {
    MAT_ASSERT(a1.cols == a2.rows);
    size_t n = a1.cols;
    MAT_ASSERT(a1.rows == dest.rows);
    MAT_ASSERT(a2.cols == dest.cols);

    for (size_t row = 0; row < dest.rows; row++) {
        for (size_t col = 0; col < dest.cols; col++) {
            MAT_AT(dest, row, col) = 0;
            for (size_t i = 0; i < n; i++) {
                MAT_AT(dest, row, col) += MAT_AT(a1, row, i) * MAT_AT(a2, i, col);
            }
        }
    }
}

Matrix __multiply__matrix(int num_args, ...) {
    va_list args;
    va_start(args, num_args);
    
    Matrix a0;
    Matrix dest;

    for (size_t i = 0; i < num_args; i++) {
        Matrix a1 = va_arg(args, Matrix);
        if(i == 0) a0 = a1;
        else {
            MAT_ASSERT(a0.cols == a1.rows);
            dest = __alloc__matrix(a0.rows, a1.cols);
            __dot__matrix(dest, a0, a1);
            a0 = dest;
        } 
    }
    
    va_end(args);
    return dest;
}

void __copy__matrix(Matrix dest, Matrix src) {
    MAT_ASSERT(dest.rows == src.rows);
    MAT_ASSERT(dest.cols == src.cols);

    for (size_t row = 0; row < dest.rows; row++) {
        for (size_t col = 0; col < dest.cols; col++) {
            MAT_AT(dest, row, col) = MAT_AT(src, row, col);
        }
    }
}

void __sigmoidf__matrix(Matrix mat) {
    for (size_t row = 0; row < mat.rows; row++) {  
        for (size_t col = 0; col < mat.cols; col++) {
            MAT_AT(mat, row, col) = sigmoidf(MAT_AT(mat, row, col));
        }
    }
}

void __relu__matrix(Matrix mat) {
    for (size_t row = 0; row < mat.rows; row++) {  
        for (size_t col = 0; col < mat.cols; col++) {
            MAT_AT(mat, row, col) = relu(MAT_AT(mat, row, col));
        }
    }
}

void __scalar__multiply(Matrix mat, float x) {
    for (size_t row = 0; row < mat.rows; row++) {
        for (size_t col = 0; col < mat.cols; col++) {
            MAT_AT(mat, row, col) *= x;
        }
    }
}

#endif