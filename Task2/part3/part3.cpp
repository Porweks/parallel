#include <iostream>
#include <omp.h>
#include <cmath>

#define N 200
#define nt 5
#define variable 0.01

void matrix_vector_product_omp(double *a, double *b, double *c, int m, int n){
    #pragma omp parallel num_threads(nt)
    {
        int nthreads = omp_get_num_threads();
        int threadid = omp_get_thread_num();
        int items_per_thread = m / nthreads;
        int lb = threadid * items_per_thread;
        int ub = (threadid == nthreads - 1) ? (m - 1) : (lb + items_per_thread - 1);
        // #pragma omp for
        for (int i = lb; i <= ub; i++) {
            c[i] = 0.0;
            for (int j = 0; j < n; j++)
                c[i] += a[i * n + j] * b[j];

        }
    }
}

double* solution_without_for(double* A, double* b){
    double *x = new double[N];
    // #pragma omp parallel num_threads(nt)
    // {
    //     int nthreads = omp_get_num_threads();
    //     int threadid = omp_get_thread_num();
    //     int items_per_thread = N / nthreads;
    //     int lb = threadid * items_per_thread;
    //     int ub = (threadid == nthreads - 1) ? (N - 1) : (lb + items_per_thread - 1);
    //     for (int i = lb; i <= ub; i++){
        for (;;){
            double *c;
            c=new double[N];
            matrix_vector_product_omp(A,x,c,N,N);
            for(int i=0; i<N;i++)
                c[i]=(c[i] - b[i]) * variable;
                
            for(int i=0; i<N;i++){
                #pragma omp atomic
                x[i]-=c[i];
            }
            double *check=new double[N];
            matrix_vector_product_omp(A,x,check,N,N);
            if()

        }
    // }
    return x;
}

double* solution_with_for(double* A, double* b){
    for(int i = 0; i < N; i++){
        double div = A[i*N+i];
        for(int j = 0; j < N; j++){
            A[i*N+j] /=div;
        }
        b[i] /= div;
        for(int k = 0; k < N; k++){
            if(k!=i){
                double mult = A[k*N+i];
                for (int l = 0; l < N; l++){
                     A[k*N+l] -= mult * A[i*N+l];
                }
                b[k] -= mult * b[i];
            }
        }
    }
    return b;
} 

int main(){
    double *A;
    A = new double[N*N];
    for(int i = 0 ; i < N ; i++)
        for(int j = 0; j < N; j++){
            A[i * N + j]=1;
            if(i == j)
                A[i * N + j]=2;
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
    printf("Elapsed time (parallel without for): %.6f sec.\n", end-start);
    for(int i=0;i<N;i++)
        printf("[%f],", x[i]);
    
    // printf("\n");
    // start = omp_get_wtime(); 
    // x = solution_with_for(A,b);
    // end = omp_get_wtime();
    // printf("Elapsed time (parallel with for): %.6f sec.\n", end-start);
    // for(int i = 0; i < N; i++)
    //     printf("[%f],", x[i]);
    // printf("\n");
    return 0;
}