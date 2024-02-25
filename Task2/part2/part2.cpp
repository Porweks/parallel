#include <iostream>
#include <omp.h>
#include <cmath>


const double a = -4.0;
const double b = 4.0;
const int nsteps = 40000000;


double func(double x)
{
    return exp(-x * x);
}


double integrate(double a, double b, int n)
{
    double h = (b - a) / n;
    double sum = 0.0;
    for (int i = 0; i < n; i++)
        sum += func(a + h * (i + 0.5));
    sum *= h;
    return sum;
}

double integrate_omp(double a, double b, int n)
{
    double h = (b - a) / n;
    double sum = 0.0;
    #pragma omp parallel num_threads(7)
    {
        int nthreads = omp_get_num_threads();
        int threadid = omp_get_thread_num();
        int items_per_thread = n / nthreads;
        int lb = threadid * items_per_thread;
        int ub = (threadid == nthreads - 1) ? (n - 1) : (lb + items_per_thread - 1);
        double sumloc = 0.0;
        for (int i = lb; i <= ub; i++)
            sumloc += func(a + h * (i + 0.5));
        #pragma omp atomic
        sum+=sumloc;
    }
    sum *= h;
    return sum;
}

int main(){

    // double start = omp_get_wtime();
    // double res_serial = integrate(a, b, nsteps);
    // double end = omp_get_wtime(); 
    // printf("Elapsed time (serial): %.6f sec.\n", end-start);
    double start = omp_get_wtime();
    double res_paralel = integrate_omp(a, b, nsteps);
    double end = omp_get_wtime(); 
    printf("Elapsed time (paralel): %.6f sec.\n", end-start);
}