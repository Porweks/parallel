#include <iostream>
#include <omp.h>
#include <cmath>

const int N = 46340;
const int nt = 80;
const double variable = (double)1/N;

void mx_vec_mult(double* mx, double* vec, double* res){
    #pragma omp for
    for (int i = 0; i < N; i++){
        int sum = 0;
        for(int j = 0; j<N; j++){
            sum += mx[i*N + j] * vec[j];
        }
        res[i] = sum;
    }
}

void vec_subtraction(double* vec_1, double* vec_2, double* res){
    #pragma omp for
    for (int i = 0; i < N; i++){
        res[i] = vec_1[i] - vec_2[i];
    }
}

void inplace_vec_subtraction(double* vec_1, double* vec_2){
    #pragma omp for
  for (int i = 0; i < N; i++){
        vec_1[i] -= vec_2[i];
    }
}

void scale_by(double* vec){
    #pragma omp for
    for (int i = 0; i < N; i++){
        vec[i] *= variable;
    }
}

double find_norm(double* vec){
    double norm = 0;
        double sumloc = 0.0;
        #pragma omp for
        for (int i = 0; i < N; i++){
            norm += pow(vec[i], 2);
        }
        #pragma omp atomic
        norm += sumloc;
    return sqrt(norm);
}

void simple_iteration(double* matrix, double* vec_x, double* vec_b){
    double* ax = (double*)malloc(N*sizeof(double));
    double* subtracted = (double*)malloc(N*sizeof(double));
    double norm_b = find_norm(vec_b);
    double norm_sub;

    #pragma omp parallel
    {
    while(1){
    
        mx_vec_mult(matrix, vec_x, ax);
        vec_subtraction(ax, vec_b, subtracted);
        norm_sub = find_norm(subtracted);
        printf("%f\n", norm_sub / norm_b);
        if (norm_sub / norm_b < 0.0000001)
            break;
        scale_by(subtracted);
        inplace_vec_subtraction(vec_x, subtracted);
    }
    }
    free(ax);
    free(subtracted);
}

int main(){
    double *A;
    A = (double*)malloc(N*N*sizeof(double));
    for(int i = 0 ; i < N ; i++)
        for(int j = 0; j < N; j++){
            if(i == j)
                A[i * N + j]=2;
            else 
                A[i * N + j]=1;
        }
    double *b;
    b = (double*)malloc(N*sizeof(double));   
    for(int i = 0; i < N; i++)
        b[i] = N + 1;
    double *x;
    x = (double*)malloc(N*sizeof(double));
    double start = omp_get_wtime(); 
    simple_iteration(A,x,b);
    double end = omp_get_wtime(); 
    free(A);
    free(x);
    free(b);
    // for(int i=0;i<N;i++)
    //     printf("[%f],", x[i]);
    printf("Elapsed time (parallel without for): %.6f sec.\n", end-start);
    return 0;
}