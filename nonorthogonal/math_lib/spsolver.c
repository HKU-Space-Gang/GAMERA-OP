#include "spsolver.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <limits.h>


/* ============ 选择默认后端：1=RCM+带状LU（默认），0=稠密LU ============ */
#ifndef SP_SOLVER_USE_BAND_RCM
#define SP_SOLVER_USE_BAND_RCM 1
#endif

/*---------------- 安全小工具 ----------------*/
static inline int safe_mul_size(size_t a, size_t b, size_t *out){
    if(a==0 || b==0){ *out=0; return 1; }
    if(a > SIZE_MAX / b) return 0;
    *out = a*b; return 1;
}

/*======================= Triplet -> CSR =======================*/
typedef struct { int i,j; double x; } TItem;
static int cmp_titem_(const void *a, const void *b){
    const TItem *A=(const TItem*)a,*B=(const TItem*)b;
    if(A->i!=B->i) return (A->i<B->i)?-1:1;
    if(A->j!=B->j) return (A->j<B->j)?-1:1;
    return 0;
}

CSR triplet_to_csr(int n, int nz, const int *Ti, const int *Tj, const double *Tx){
    CSR A = (CSR){0}; A.n = n;
    TItem *arr = (TItem*)malloc(sizeof(TItem)*nz);
    if(!arr){ fprintf(stderr,"triplet_to_csr: OOM\n"); exit(1); }
    for(int k=0;k<nz;++k){ arr[k].i=Ti[k]; arr[k].j=Tj[k]; arr[k].x=Tx[k]; }
    qsort(arr, nz, sizeof(TItem), cmp_titem_);

    int *row_nnz = (int*)calloc(n, sizeof(int));
    if(!row_nnz){ fprintf(stderr,"triplet_to_csr: OOM\n"); exit(1); }

    int m = 0;
    for(int k=0;k<nz;){
        int i=arr[k].i, j=arr[k].j; double x=arr[k].x;
        int k2=k+1;
        while(k2<nz && arr[k2].i==i && arr[k2].j==j){ x+=arr[k2].x; ++k2; }
        if(x!=0.0){ ++row_nnz[i]; ++m; }
        k=k2;
    }

    A.nnz = m;
    A.Ap = (int*)malloc(sizeof(int)*(n+1));
    A.Aj = (int*)malloc(sizeof(int)*m);
    A.Ax = (double*)malloc(sizeof(double)*m);
    if(!A.Ap||!A.Aj||!A.Ax){ fprintf(stderr,"triplet_to_csr: OOM\n"); exit(1); }

    A.Ap[0]=0; for(int i=0;i<n;++i) A.Ap[i+1]=A.Ap[i]+row_nnz[i];
    int *cursor = (int*)malloc(sizeof(int)*n);
    if(!cursor){ fprintf(stderr,"triplet_to_csr: OOM\n"); exit(1); }
    memcpy(cursor, A.Ap, sizeof(int)*n);

    for(int k=0;k<nz;){
        int i=arr[k].i, j=arr[k].j; double x=arr[k].x;
        int k2=k+1;
        while(k2<nz && arr[k2].i==i && arr[k2].j==j){ x+=arr[k2].x; ++k2; }
        if(x!=0.0){
            int p = cursor[i]++;
            A.Aj[p] = j; A.Ax[p] = x;
        }
        k=k2;
    }
    free(arr); free(row_nnz); free(cursor);
    return A;
}

void free_csr(CSR *A){
    if(!A) return;
    free(A->Ap); free(A->Aj); free(A->Ax);
    A->Ap = A->Aj = NULL; A->Ax = NULL;
    A->n = A->nnz = 0;
}

/*======================= CSR -> 稠密 =======================*/
static void csr_to_dense_rowmajor(const CSR *A, double *M){ // M: n*n
    int n=A->n;
    for(int i=0;i<n*n;++i) M[i]=0.0;
    for(int i=0;i<n;++i){
        for(int p=A->Ap[i]; p<A->Ap[i+1]; ++p){
            int j=A->Aj[p];
            M[i*(size_t)n + j] = A->Ax[p];
        }
    }
}

/*======================= 稠密 LU(部分选主元) =======================*/
static int lu_factor_pp(double *LU, int n, int *piv){
    for(int i=0;i<n;++i) piv[i]=i;
    for(int k=0;k<n;++k){
        int pivrow=k; double amax=fabs(LU[k*(size_t)n+k]);
        for(int r=k+1;r<n;++r){
            double v=fabs(LU[r*(size_t)n+k]);
            if(v>amax){ amax=v; pivrow=r; }
        }
        if(amax==0.0) return -1; /* 奇异 */

        if(pivrow!=k){
            for(int c=0;c<n;++c){
                double t=LU[k*(size_t)n+c];
                LU[k*(size_t)n+c]=LU[pivrow*(size_t)n+c];
                LU[pivrow*(size_t)n+c]=t;
            }
            int tp=piv[k]; piv[k]=pivrow; piv[pivrow]=tp;
        }

        for(int i=k+1;i<n;++i){
            LU[i*(size_t)n+k] /= LU[k*(size_t)n+k];
            double lik=LU[i*(size_t)n+k];
            for(int j=k+1;j<n;++j) LU[i*(size_t)n+j] -= lik*LU[k*(size_t)n+j];
        }
    }
    return 0;
}

static void lu_solve_pp(const double *LU, const int *piv, const double *b, double *x, int n){
    for(int i=0;i<n;++i) x[i]=b[piv[i]];                /* Pb */
    for(int i=0;i<n;++i) for(int j=0;j<i;++j) x[i]-=LU[i*(size_t)n+j]*x[j]; /* Ly */
    for(int i=n-1;i>=0;--i){
        for(int j=i+1;j<n;++j) x[i]-=LU[i*(size_t)n+j]*x[j];
        x[i] /= LU[i*(size_t)n+i];
    }
}

int solve_from_csr_dense_lu(const CSR *A, const double *b, double *x){
    int n=A->n;
    size_t nelem;
    if(!safe_mul_size((size_t)n,(size_t)n,&nelem)) return -2;
    double *M=(double*)malloc(sizeof(double)*nelem);
    int *piv=(int*)malloc(sizeof(int)*n);
    if(!M||!piv){ free(M); free(piv); return -2; }
    csr_to_dense_rowmajor(A,M);
    if(lu_factor_pp(M,n,piv)!=0){ free(M); free(piv); return -1; }
    lu_solve_pp(M,piv,b,x,n);
    free(M); free(piv);
    return 0;
}

/*======================= RCM 重排序（基于无向图） =======================*/
typedef struct { int *d; int sz, cap; } IntVec;
static void iv_init(IntVec *v){ v->d=NULL; v->sz=0; v->cap=0; }
static void iv_push(IntVec *v,int x){
    if(v->sz==v->cap){
        v->cap = v->cap? v->cap*2:8;
        v->d = (int*)realloc(v->d,sizeof(int)*v->cap);
    }
    v->d[v->sz++]=x;
}
static int cmp_int_(const void*a,const void*b){ int x=*(const int*)a,y=*(const int*)b; return (x>y)-(x<y); }

static void build_undirected_adj_from_csr(const CSR*A, IntVec **adj, int *deg){
    int n=A->n;
    *adj=(IntVec*)malloc(sizeof(IntVec)*n);
    for(int i=0;i<n;++i){ iv_init(&(*adj)[i]); deg[i]=0; }
    for(int i=0;i<n;++i){
        for(int p=A->Ap[i]; p<A->Ap[i+1]; ++p){
            int j=A->Aj[p];
            if(i==j) continue;
            iv_push(&(*adj)[i], j);
            iv_push(&(*adj)[j], i);
        }
    }
    for(int i=0;i<n;++i){
        if((*adj)[i].sz==0){ deg[i]=0; continue; }
        qsort((*adj)[i].d, (*adj)[i].sz, sizeof(int), cmp_int_);
        int w=0;
        for(int r=0;r<(*adj)[i].sz;++r){
            if(w==0 || (*adj)[i].d[r]!=(*adj)[i].d[w-1]) (*adj)[i].d[w++]=(*adj)[i].d[r];
        }
        (*adj)[i].sz=w; deg[i]=w;
    }
}

static void rcm_order(const CSR*A, int **perm_out, int **invperm_out){
    int n=A->n;
    IntVec *adj=(IntVec*)malloc(sizeof(IntVec)*n);
    int *deg=(int*)malloc(sizeof(int)*n);
    build_undirected_adj_from_csr(A,&adj,deg);

    int *visited=(int*)calloc(n,sizeof(int));
    int *Q=(int*)malloc(sizeof(int)*n); int qh=0, qt=0;
    int *R=(int*)malloc(sizeof(int)*n); int rsz=0;

    while(rsz<n){
        int s=-1, mind=INT_MAX;
        for(int i=0;i<n;++i) if(!visited[i] && deg[i]<mind){ mind=deg[i]; s=i; }
        visited[s]=1; Q[qt++]=s;
        while(qh<qt){
            int u=Q[qh++]; R[rsz++]=u;
            int sz=adj[u].sz;
            if(sz){
                int *ord=(int*)malloc(sizeof(int)*sz);
                memcpy(ord, adj[u].d, sizeof(int)*sz);
                for(int i=1;i<sz;++i){
                    int key=ord[i], j=i-1;
                    while(j>=0 && deg[ord[j]]>deg[key]){ ord[j+1]=ord[j]; --j; }
                    ord[j+1]=key;
                }
                for(int i=0;i<sz;++i){
                    int v=ord[i];
                    if(!visited[v]){ visited[v]=1; Q[qt++]=v; }
                }
                free(ord);
            }
        }
    }

    int *perm=(int*)malloc(sizeof(int)*n), *inv=(int*)malloc(sizeof(int)*n);
    for(int i=0;i<n;++i){ int old=R[n-1-i]; perm[i]=old; }
    for(int newi=0; newi<n; ++newi) inv[ perm[newi] ] = newi;

    for(int i=0;i<n;++i) free(adj[i].d);
    free(adj); free(deg); free(visited); free(Q); free(R);

    *perm_out=perm; *invperm_out=inv;
}

/*======================= 对称置换 B = P*A*P^T =======================*/
static CSR csr_sym_permute(const CSR*A, const int *invperm){
    int n=A->n, nz=A->nnz, k=0;
    int *Ti=(int*)malloc(sizeof(int)*nz);
    int *Tj=(int*)malloc(sizeof(int)*nz);
    double *Tx=(double*)malloc(sizeof(double)*nz);
    for(int i=0;i<n;++i){
        for(int p=A->Ap[i]; p<A->Ap[i+1]; ++p){
            int j=A->Aj[p]; double x=A->Ax[p];
            Ti[k]=invperm[i]; Tj[k]=invperm[j]; Tx[k]=x; ++k;
        }
    }
    CSR B=triplet_to_csr(n,nz,Ti,Tj,Tx);
    free(Ti); free(Tj); free(Tx);
    return B;
}

/*======================= 计算带宽 =======================*/
static int csr_bandwidth(const CSR*A){
    int n=A->n, bw=0;
    for(int i=0;i<n;++i){
        for(int p=A->Ap[i]; p<A->Ap[i+1]; ++p){
            int j=A->Aj[p];
            int d = (i>j)? (i-j):(j-i);
            if(d>bw) bw=d;
        }
    }
    return bw;
}

/*======================= CSR -> 带状存储 =======================*/
/* 带状存储：AB 行数 = 2*bw+1, 列数 = n
   下标: AB[ bw + (i-j) ][ j ] = A(i,j), 仅当 |i-j|<=bw */
static inline double* AB_at(double *AB, int n, int bw, int i, int j){
    int r = bw + (i - j);
    if(r<0 || r>2*bw) return NULL;
    return &AB[ ((size_t)r)*(size_t)n + (size_t)j ];
}

static double* csr_to_band(const CSR *A, int bw){
    if(bw<0){ fprintf(stderr,"csr_to_band: negative bw\n"); return NULL; }
    size_t rows = (size_t)bw*2u + 1u;
    size_t cols = (size_t)A->n;
    size_t nelem;

    if (!safe_mul_size((size_t)rows, (size_t)cols, &nelem)) {
        fprintf(stderr, "[csr_to_band] band matrix too large, overflow in ldab*n\n");
        return NULL;    // 出错时返回 NULL，而不是 -1
    }

    double *AB = (double*)calloc(nelem, sizeof(double));
    
    if (!AB) {
        fprintf(stderr, "[csr_to_band] calloc failed for band matrix (nelem=%zu)\n", nelem);
        return NULL;    // 同样返回 NULL
    }

    for(int i=0;i<A->n;++i){
        for(int p=A->Ap[i]; p<A->Ap[i+1]; ++p){
            int j=A->Aj[p];
            int r=bw+(i-j);
            if(r<0 || r>2*bw) continue;
            AB[ ((size_t)r)*cols + (size_t)j ] += A->Ax[p];
        }
    }
    return AB;
}

/*======================= 带状 LU（无主元） =======================*/
static int band_lu_nopiv(double *AB, int n, int bw){
    const double eps = 1e-30;
    for(int k=0;k<n;++k){
        double *Ukkp = AB_at(AB,n,bw,k,k);
        if(!Ukkp) return -2;
        double Ukk = *Ukkp;
        if(fabs(Ukk) < eps) return -1;

        int i_max = (k+bw < n)? (k+bw) : (n-1);
        int j_max = i_max;

        for(int i=k+1; i<=i_max; ++i){
            double *Likp = AB_at(AB,n,bw,i,k);
            if(Likp) *Likp /= Ukk;
        }
        for(int i=k+1; i<=i_max; ++i){
            double *Likp = AB_at(AB,n,bw,i,k);
            if(!Likp) continue;
            double Lik = *Likp;
            if(Lik==0.0) continue;
            for(int j=k+1; j<=j_max; ++j){
                double *Ukjp = AB_at(AB,n,bw,k,j);
                if(!Ukjp) continue;
                double *Aijp = AB_at(AB,n,bw,i,j);
                if(Aijp) *Aijp -= Lik * (*Ukjp);
            }
        }
    }
    return 0;
}

/*======================= 带状三角求解 =======================*/
static void band_forward_subst(double *AB, int n, int bw, double *b){ // L unit diag
    for(int i=0;i<n;++i){
        int j0 = (i-bw>0)? (i-bw):0;
        for(int j=j0; j<i; ++j){
            double *Lij = AB_at(AB,n,bw,i,j);
            if(Lij) b[i] -= (*Lij)*b[j];
        }
    }
}
static void band_backward_subst(double *AB, int n, int bw, double *x){ // U diag non-unit
    for(int i=n-1;i>=0; --i){
        int j1 = (i+bw < n)? (i+bw) : (n-1);
        for(int j=i+1; j<=j1; ++j){
            double *Uij = AB_at(AB,n,bw,i,j);
            if(Uij) x[i] -= (*Uij)*x[j];
        }
        double *Uii = AB_at(AB,n,bw,i,i);
        x[i] /= *Uii;
    }
}

/*======================= 顶层：RCM + 带状LU求解 =======================*/
void destroy_lu_band(LU_Band *F){
    if(!F) return; free(F->perm); free(F->invperm); free(F->AB);
    F->perm=F->invperm=NULL; F->AB=NULL; F->n=F->bw=0;
}

int factorize_lu_band_rcm(const CSR*A, LU_Band *F){
    memset(F,0,sizeof(*F));
    int *perm=NULL,*invperm=NULL;
    rcm_order(A,&perm,&invperm);
    CSR B = csr_sym_permute(A, invperm);
    int bw = csr_bandwidth(&B);
    double *AB = csr_to_band(&B, bw);
    if(!AB){ free(perm); free(invperm); free_csr(&B); return -2; }
    int rc = (AB? band_lu_nopiv(AB, B.n, bw) : -2);
    free_csr(&B);
    if(rc!=0){ free(perm); free(invperm); free(AB); return rc; }
    F->n = A->n; F->perm=perm; F->invperm=invperm; F->bw=bw; F->AB=AB;
    return 0;
}

int solve_lu_band_rcm(const LU_Band *F, const double *b, double *x){
    int n=F->n, bw=F->bw;
    double *bp=(double*)malloc(sizeof(double)*n);
    if(!bp) return -2;
    /* b' = P*b  (perm[new]=old) */
    for(int newi=0; newi<n; ++newi){ int old=F->perm[newi]; bp[newi]=b[old]; }
    band_forward_subst(F->AB, n, bw, bp);
    band_backward_subst(F->AB, n, bw, bp);
    /* x = P^T x' */
    for(int newi=0; newi<n; ++newi){ int old=F->perm[newi]; x[old]=bp[newi]; }
    free(bp);
    return 0;
}

/*======================= 包装：solve_linear_system =======================*/
int solve_linear_system(CSR *A, double *b, double *x){
#if SP_SOLVER_USE_BAND_RCM
    LU_Band F; int rc = factorize_lu_band_rcm(A, &F);
    if(rc!=0) return rc;
    rc = solve_lu_band_rcm(&F, b, x);
    destroy_lu_band(&F);
    return rc;
#else
    return solve_from_csr_dense_lu(A, b, x);
#endif
}
