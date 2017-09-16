#include <iostream>
#include <cmath>
//#include "jacobi.h"



using namespace std;

void printMatrixA(double** A, int N);

int main(){
    //Definerer var som trengs for a lage matrisen A
    int N = 5;
    double rho_0 = 0;
    double rho_max = 5;
    double h = (rho_max - rho_0)/((double)N);
    double hh = h*h;
    double e = -1/hh;

    //Definerer var som trengs for Jacobi-metoden:
    int k, l;
    double epsilon = 1.0e-8;
    double max_number_iterations = (double) N * (double) N * (double) N;
    int iterations = 0;
    //double max_offdiag = maxoffdiag(A, &k, &l, N);


    //Lager vektorer som trengs:
    double* rho = new double[N];
    double* V = new double[N];
    double* d = new double[N];

   for (int i = 0; i < N; i++){
        rho[i] = rho_0 + i*h;
        V[i] = rho[i]*rho[i];
        d[i] = 2./hh + V[i];
   }

    //Lager matrisen A:
    double** A = new double*[N];
    for (int i = 0; i < N; i++){
        A[i] = new double[N];
    }

    for (int i = 0; i < N; i++){
        for (int j = 0; j < N; j++){
            A[i][i] = d[i];
            if (i == j-1){
                A[i][j] = e;
            }
            else if (i == j+1){
                    A[i][j] = e;
            }
            else{
                A[i][j] = 0.0;
            }
        }
    }

    printMatrixA(A, N);

    return 0;
}

void printMatrixA(double** A, int N) {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            printf("%6.2f ", A[i][j]);
        }
        cout << endl;
    }
    cout << endl;
    fflush(stdout);
}

//Setting up the eigenvector matrix:
void jacobi_method(double** A, double** R, int N){
    for (int i = 0; i < N; i++){
        for (int j = 0; j < N; j++){
            if (i == j){
                R[i][j] = 1.0;
            }
            else{
                R[i][j] = 0.0;
            }
        }
    }
}

//Function to find the maximum element:
double maxoffdiag ( double ** A, int * k, int * l, int n ){
    double max = 0.0;
    for ( int i = 0; i < n; i++ ) {
        for ( int j = i + 1; j < n; j++ ) {
            if ( fabs(A[i][j]) > max ){
                max = fabs(A[i][j]);
                *l = i;
                *k = j;
            }
        }
    }
return max;
}

