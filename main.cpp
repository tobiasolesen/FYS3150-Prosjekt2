#include <iostream>

using namespace std;

void printMatrixA(double** A, int N);

int main(){
    int N = 5;
    double rho_0 = 0;
    double rho_max = 5;
    double h = (rho_max - rho_0)/((double)N);
    double hh = h*h;
    double e = -1/hh;

    //Lager vektoren rho:
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
