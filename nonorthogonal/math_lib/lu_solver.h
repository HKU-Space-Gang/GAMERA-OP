/**
 * lu_solver.h
 * 
 * LU分解求解线性方程组的声明头文件
 * LU decomposition and solver (declarations only)
 */

#ifndef LU_SOLVER_H
#define LU_SOLVER_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/**
 * LU分解（带部分主元选取）LU decomposition with partial pivoting
 * 
 * @param A    n×n矩阵（输入/输出，分解后L和U存储在A中） n×n matrix (input/output, L and U stored in A after decomposition)
 * @param perm 排列数组（记录行交换） permutation array (records row exchanges)
 * @param n    矩阵维度 matrix dimension
 * @return     0成功，-1失败（矩阵奇异） 0 on success, -1 on failure (singular matrix)
 */
int lu_decomp(double **A, int *perm, int n);

/**
 * 求解 Ax = b，其中A已经LU分解 Solve Ax = b where A is already LU decomposed
 * 
 * @param A    LU分解后的矩阵 LU decomposed matrix
 * @param perm 排列数组 permutation array
 * @param b    右端向量 right-hand side vector
 * @param x    解向量（输出） solution vector (output)
 * @param n    维度 dimension
 */
void lu_solve(double **A, int *perm, double *b, double *x, int n);

/**
 * 分配n×n矩阵 Allocate n×n matrix
 * 
 * @param n 矩阵维度 matrix dimension
 * @return  指向矩阵的指针  pointer to the matrix
 */
double** lu_alloc_matrix(int n);

/**
 * 释放矩阵内存 Free matrix memory
 * 
 * @param mat 矩阵指针 matrix pointer
 * @param n   矩阵维度 matrix dimension
 */
void lu_free_matrix(double **mat, int n);

/**
 * 复制矩阵（因为LU分解会修改原矩阵） Copy matrix (because LU decomposition modifies the original matrix)
 * 
 * @param src 源矩阵 source matrix
 * @param dst 目标矩阵 destination matrix
 * @param n   矩阵维度   matrix dimension
 */
void lu_copy_matrix(double **src, double **dst, int n);

/**
 * 分配向量 Allocate vector
 * 
 * @param n 向量维度 vector dimension
 * @return  指向向量的指针 pointer to the vector
 */
double* lu_alloc_vector(int n);

/**
 * 释放向量内存 Free vector memory
 * 
 * @param vec 向量指针 vector pointer
 */
void lu_free_vector(double *vec);

/**
 * 分配排列数组 Allocate permutation array
 * 
 * @param n 数组大小 array size
 * @return  指向排列数组的指针 pointer to the permutation array
 */
int* lu_alloc_permutation(int n);

/**
 * 释放排列数组内存 Free permutation array memory
 * 
 * @param perm 排列数组指针 permutation array pointer
 */
void lu_free_permutation(int *perm);

/**
 * 打印矩阵（用于调试） Print matrix (for debugging)
 * 
 * @param A 矩阵指针 matrix pointer
 * @param n 矩阵维度 matrix dimension
 */
void lu_print_matrix(double **A, int n);

/**
 * 打印向量（用于调试） Print vector (for debugging)
 * 
 * @param v 向量指针 vector pointer
 * @param n 向量维度 vector dimension
 */
void lu_print_vector(double *v, int n);

#endif /* LU_SOLVER_H */
