#ifndef SP_SOLVER_H_
#define SP_SOLVER_H_

#ifdef __cplusplus
extern "C" {
#endif

/* CSR */
typedef struct {
    int n, nnz;
    int *Ap, *Aj;
    double *Ax;
} CSR;

CSR  triplet_to_csr(int n, int nz, const int *Ti, const int *Tj, const double *Tx);
void free_csr(CSR *A);

/* 先声明 LU_Band！ */
typedef struct {
    int n;
    int *perm, *invperm;
    int bw;
    double *AB;
} LU_Band;

/* 稠密LU */
int solve_from_csr_dense_lu(const CSR *A, const double *b, double *x);

/* 带状LU + RCM */
int  factorize_lu_band_rcm(const CSR *A, LU_Band *F);
int  solve_lu_band_rcm(const LU_Band *F, const double *b, double *x);
void destroy_lu_band(LU_Band *F);

/* 统一入口（可选） */
int solve_linear_system(CSR *A, double *b, double *x);

#ifdef __cplusplus
}
#endif
#endif /* SP_SOLVER_H_ */
