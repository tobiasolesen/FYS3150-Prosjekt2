//Bruk &k nar du vil endre variabelen k inni funksjonen. Brukes i prototypen av funksjonen og i def av funksjonen.
//Ikke i kallet pa funksjonen. Trenger da ikke bruke *k noe sted tror jeg.
//Matriser skal vaere dobbeltpekere.
//Vektorer skal vaere enkeltpekere.
//Vet at de 3 laveste egenverdiene skal vaere: 3, 7, 11

#include <iostream>
#include <cmath>
#include <array>
#include <vector>
//#include "jacobi.h"
using namespace std;

//Function prototypes:
void printMatrixA(double** A, int N);
double maxoffdiag(double** A, int& k, int& l, int N);
void jacobi_method(double** A, double** R, int N);
void rotate ( double ** A, double ** R, int k, int l, int N );

int main(){
    //Definerer var som trengs for a lage matrisen A
    int N = 200;
    double rho_0 = 0;
    double rho_max = 5;
    double h = (rho_max - rho_0)/((double)N);
    double hh = h*h;
    double e = -1/hh;

    //Lager vektorer som trengs:
    double* rho = new double[N];
    double* V = new double[N];
    double* d = new double[N];

    //Fyller vektorene
    for (int i = 0; i < N; i++){
        rho[i] = rho_0 + (i+1)*h;
        V[i] = rho[i]*rho[i];
        d[i] = 2./hh + V[i];
    }

    //Deklarerer/lager matrisene A og R:
    double** A = new double*[N];
    double** R = new double*[N];
    for (int i = 0; i < N; i++){
        A[i] = new double[N];
        R[i] = new double[N];
    }

    //Fyller inn offdiag-elementene i A:
    for (int i = 0; i < N; i++){
        for (int j = 0; j < N; j++){
            if (i == j-1){
                A[i][j] = e;
            } else if (i == j+1){
                A[i][j] = e;
            } else {
                A[i][j] = 0.0;
            }
        }
    }
    //Fyller diagonalelementene i A:
    for (int i = 0; i < N; i++){
        A[i][i] = d[i];
    }

    //printMatrixA(A, N);  //Skriver ut matrisen A som vi starter med

    jacobi_method(A, R, N);

    return 0;
}  //Slutt pa main

//Funksjon som skriver ut matrise
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

//Setting up the eigenvector matrix (som identitetsmatrisen?):
void jacobi_method(double** A, double** R, int N) {

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

    int k, l;
    double epsilon = 1.0e-8;
    double max_number_iterations = (double) N * (double) N * (double) N;
    int iterations = 0;
    double max_offdiag = maxoffdiag(A, k, l, N);

    while ( fabs(max_offdiag) > epsilon && (double) iterations < max_number_iterations ){
        max_offdiag = maxoffdiag ( A, k, l, N );
        rotate ( A, R, k, l, N );
        iterations++;
        //Hopper ut av denne lokken nar offdiag elementene til A er essentially zero
    }

    cout << "Number of iterations: " << iterations << "\n";

    vector<double>diag(N);  //Deklarerer vektor diag med lengde N?

    for (int i = 0; i < N; i++){
        diag[i] = A[i][i];  //Fyller diag med A sine diagonalelementer
    }

    sort(diag.begin(), diag.end());  //Sorterer vektoren fra minst til storste verdier
    //Skriver ut diag sine 3 forste elementer (3 forste egenverdiene):
    cout << "3 forste egenverdier: " << diag[0] << " " << diag[1] << " " << diag[2] << endl;

    return;
}

//Function to find the maximum element of A:
double maxoffdiag ( double** A, int& k, int& l, int N) {

    double max = 0.0;
    for ( int i = 0; i < N; i++ ) {
        for ( int j = i + 1; j < N; j++ ) {
            if ( fabs(A[i][j]) > max ){
                max = fabs(A[i][j]);
                l = i;
                k = j;  //Maks-elementet girs indeksene [k][l], og [l][k] senere?
            }
        }
    }

    return max;

}

// Function to find the values of cos and sin
void rotate ( double ** A, double ** R, int k, int l, int N ){

    double s, c;
    //cout << k << " " << l << endl;

    if ( A[k][l] != 0.0 ) {

        double t, tau;
        tau = (A[l][l] - A[k][k])/(2*A[k][l]);

        if ( tau > 0 ) {
            t = 1.0/(tau + sqrt(1.0 + tau*tau));
        }
        else {
            t = -1.0/( -tau + sqrt(1.0 + tau*tau));
        }
        c = 1.0/sqrt(1.0+t*t);
        s = c*t;
    } else{
        c = 1.0;
        s = 0.0;
    }

    double a_kk, a_ll, a_ik, a_il, r_ik, r_il;
    a_kk = A[k][k];
    a_ll = A[l][l];

    // changing the matrix elements with indices k and l A[k][k] = c*c*a_kk - 2.0*c*s*A[k][l] + s*s*a_ll; A[l][l] = s*s*a_kk + 2.0*c*s*A[k][l] + c*c*a_ll; A[k][l] = 0.0; // hard-coding of the zeros

    A[k][k] = a_kk*c*c - 2*A[k][l]*c*s + a_ll*s*s;
    A[l][l] = a_ll*c*c + 2*A[k][l]*c*s + a_kk*s*s;
    A[l][k] = 0.0;
    A[k][l] = 0.0;
    // and then we change the remaining elements
    for ( int i = 0; i < N; i++ ) {
        if ( i != k && i != l ){
            a_ik = A[i][k];
            a_il = A[i][l];
            A[i][k] = c*a_ik - s*a_il;
            A[k][i] = A[i][k];
            A[i][l] = c*a_il + s*a_ik;
            A[l][i] = A[i][l];
        }
        // Finally, we compute the new eigenvectors
        r_ik = R[i][k];
        r_il = R[i][l];
        R[i][k] = c*r_ik - s*r_il;
        R[i][l] = c*r_il + s*r_ik;
    }

    return;

}
