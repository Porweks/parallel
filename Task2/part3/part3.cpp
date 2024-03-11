#include <iostream>
#include <omp.h>
#include <cmath>

#define N 46340
#define nt 80
#define variable (double)1/N

void matrix_vector_product_omp(double *a, double *b, double *c, int m, int n){
    #pragma omp parallel
    {
        int nthreads = omp_get_num_threads();
        int threadid = omp_get_thread_num();
        int items_per_thread = m / nthreads;
        int lb = threadid * items_per_thread;
        int ub = (threadid == nthreads - 1) ? (m - 1) : (lb + items_per_thread - 1);
        for (int i = lb; i <= ub; i++) {
            c[i] = 0.0;
            for (int j = 0; j < n; j++)
                c[i] += a[i * n + j] * b[j];

        }
    }
}

double find_norm(double* vec){
    double norm = 0;
        double sumloc = 0.0;
        #pragma omp for schedule(static)
        for (int i = 0; i < N; i++){
            norm += pow(vec[i], 2);
        }
        #pragma omp atomic
        norm += sumloc;
    return sqrt(norm);
}

double* solution_without_for(double* A, double* b){
    double *x = new double[N];

        double norm_b = find_norm(b);
        double norm_sub;
        double *sub;
        sub=new double[N];
        while(true){

            double *c;
            c=new double[N];

            matrix_vector_product_omp(A,x,c,N,N);

            #pragma omp for schedule(static)
            for(int i=0; i<N;i++){
                sub[i]=c[i] - b[i];
            }

            norm_sub=find_norm(sub);

            // std::cout<<norm_sub<<'/'<<norm_b<<'='<<norm_sub/norm_b<<'\n';
            // std::cout<<norm_sub<<'/'<<norm_b<<'='<<norm_sub/norm_b<<'\n';


            if(norm_sub/norm_b<0.0000001)
                break;
            
            #pragma omp for schedule(static)
            for (int i = 0; i < N; i++)
                sub[i] *= variable;

            #pragma omp for schedule(static)
            for (int i = 0; i < N; i++)
                x[i] -= sub[i];
                
        }
    // }
    return x;
}

int main(){
    double *A;
    A = new double[N*N];
    for(int i = 0 ; i < N ; i++)
        for(int j = 0; j < N; j++){
            if(i == j)
                A[i * N + j]=2;
            else 
                A[i * N + j]=1;
        }
    double *b;
    b = new double[N];   
    for(int i = 0; i < N; i++)
        b[i] = N + 1;
    double *x;
    x = new double[N];
    double start = omp_get_wtime(); 
    x = solution_without_for(A,b);
    double end = omp_get_wtime(); 
    // for(int i=0;i<N;i++)
    //     printf("[%f],", x[i]);
    printf("Elapsed time (parallel without for): %.6f sec.\n", end-start);
    return 0;
}