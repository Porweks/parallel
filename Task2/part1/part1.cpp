#include <iostream>
#include <omp.h>


#define size 20000

void matrix_vector_product(double *a, double *b, double *c, int m, int n)
{
    for (int i = 0; i < m; i++) {
        c[i] = 0.0;
        for (int j = 0; j < n; j++)
            c[i] += a[i * n + j] * b[j];

    }
}

void matrix_vector_product_omp(double *a, double *b, double *c, int m, int n){
    #pragma omp parallel num_threads(7)
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

int main(){
    double *a, *b, *c;
// Allocate memory for 2-d array a[m, n]
    int n = size;
    int m = size;
    a = new double[n*m];
    b = new double[n];
    c = new double[m];
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++)
            a[i * n + j] = i + j;
    }
    for (int j = 0; j < n; j++)
        b[j] = j;
    double start = omp_get_wtime(); 
    matrix_vector_product_omp(a, b, c, m, n);
    double end = omp_get_wtime(); 
    printf("Elapsed time (parallel): %.6f sec.\n", end-start);
    // start = omp_get_wtime(); 
    // matrix_vector_product(a, b, c, m, n);
    // end = omp_get_wtime(); 
    // printf("Elapsed time (parallel): %.6f sec.\n", end-start);
    free(a);
    free(b);
    free(c);
}