/**
 * lu_solver.c
 * 
 * LU分解求解线性方程组的实现文件
 * LU decomposition and solver implementation
 */

#include "lu_solver.h"

int lu_decomp(double **A, int *perm, int n) {
    int i, j, k, max_row;
    double max_val, tmp, factor;
    
    // 初始化排列数组 Initialize permutation array
    for (i = 0; i < n; i++) {
        perm[i] = i;
    }
    
    // LU分解主循环 LU decomposition main loop
    for (k = 0; k < n - 1; k++) {
        // 部分主元选取 Partial pivoting
        max_val = fabs(A[k][k]);
        max_row = k;
        for (i = k + 1; i < n; i++) {
            if (fabs(A[i][k]) > max_val) {
                max_val = fabs(A[i][k]);
                max_row = i;
            }
        }
        
        // 检查奇异性 Check for singularity
        if (max_val < 1e-12) {
            return -1;  // 矩阵奇异 Singular matrix
        }
        
        // 交换行 Swap rows
        if (max_row != k) {
            for (j = 0; j < n; j++) {
                tmp = A[k][j];
                A[k][j] = A[max_row][j];
                A[max_row][j] = tmp;
            }
            // 更新排列 Update permutation
            int tmp_perm = perm[k];
            perm[k] = perm[max_row];
            perm[max_row] = tmp_perm;
        }
        
        // 计算乘数并更新子矩阵 Compute multipliers and update submatrix
        for (i = k + 1; i < n; i++) {
            factor = A[i][k] / A[k][k];
            A[i][k] = factor;  // 存储L的元素 Store L element
            
            for (j = k + 1; j < n; j++) {
                A[i][j] -= factor * A[k][j];
            }
        }
    }
    
    // 检查最后一个对角元素 Check last diagonal element
    if (fabs(A[n-1][n-1]) < 1e-12) {
        return -1;
    }
    
    return 0;
}

void lu_solve(double **A, int *perm, double *b, double *x, int n) {
    int i, j;
    double *y = (double *)malloc(n * sizeof(double));
    double *pb = (double *)malloc(n * sizeof(double));
    
    if (!y || !pb) {
        fprintf(stderr, "lu_solve: memory allocation failed\n");
        free(y);
        free(pb);
        return;
    }
    
    // 应用排列到b，得到Pb. Apply permutation to b to get Pb
    for (i = 0; i < n; i++) {
        pb[i] = b[perm[i]];
    }
    
    // 前向替换求解 Ly = Pb. Forward substitution to solve Ly = Pb
    for (i = 0; i < n; i++) {
        y[i] = pb[i];
        for (j = 0; j < i; j++) {
            y[i] -= A[i][j] * y[j];
        }
    }
    
    // 后向替换求解 Ux = y. Backward substitution to solve Ux = y
    for (i = n - 1; i >= 0; i--) {
        x[i] = y[i];
        for (j = i + 1; j < n; j++) {
            x[i] -= A[i][j] * x[j];
        }
        x[i] /= A[i][i];
    }
    
    free(y);
    free(pb);
}

double** lu_alloc_matrix(int n) {
    double **mat = (double **)malloc(n * sizeof(double *));
    if (!mat) return NULL;

    for (int i = 0; i < n; i++) {
        mat[i] = (double *)malloc(n * sizeof(double));
        if (!mat[i]) {
            // 回滚已分配部分 Roll back partial allocations
            for (int k = 0; k < i; k++) {
                free(mat[k]);
            }
            free(mat);
            return NULL;
        }
    }
    return mat;
}

void lu_free_matrix(double **mat, int n) {
    if (!mat) return;
    for (int i = 0; i < n; i++) {
        free(mat[i]);
    }
    free(mat);
}

void lu_copy_matrix(double **src, double **dst, int n) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            dst[i][j] = src[i][j];
        }
    }
}

double* lu_alloc_vector(int n) {
    return (double *)malloc(n * sizeof(double));
}

void lu_free_vector(double *vec) {
    free(vec);
}

int* lu_alloc_permutation(int n) {
    return (int *)malloc(n * sizeof(int));
}

void lu_free_permutation(int *perm) {
    free(perm);
}

void lu_print_matrix(double **A, int n) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            printf("%10.4f ", A[i][j]);
        }
        printf("\n");
    }
}

void lu_print_vector(double *v, int n) {
    for (int i = 0; i < n; i++) {
        printf("%10.4f ", v[i]);
    }
    printf("\n");
}
