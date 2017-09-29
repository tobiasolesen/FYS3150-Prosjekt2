//Bruk &k nar du vil endre variabelen k inni funksjonen. Brukes i prototypen av funksjonen og i def av funksjonen.
//Ikke i kallet pa funksjonen. Trenger da ikke bruke *k noe sted tror jeg.
//Matriser skal vaere dobbeltpekere.
//Vektorer skal vaere enkeltpekere.
//Vet at de 3 laveste egenverdiene skal vaere: 3, 7, 11 (i b) )

//i d) skal jeg bruke koden fra b), endrer bare pa potensialet V(rho)
//rho er avstanden mellom partiklene (?)
//Vil finne den minste egenverdien (tilsvarer energien til ground state) og tilhorende egenvektor (bolgefunksjonen i ground state)

//Sammenligner energien vi faar for grunntilstand ( energi=1.25) med n=2-tilfellet i artikkelen, altsa energi=0.625 (hvor de har w_r=0.25).
//Vi far dobbelt saa mye energi (1.25 = 2*0.625).
//Dette er pa grunn av skaleringen av Schrodinger ligningen. Siden de beholder faktoren 1/2 i SL sa er det riktig at jeg faar dobbelt saa mye.

//apne tekstfil i terminal: less "filnavn"
//Komme ut av fil og tilbake i terminal: shift + q
//Med w_r=0.25 skal de minste egenverdiene bli: 1.25, 2.19, 3.18 ca
//Normen er bevart under similarity transform
// 10.986 rundes opp sa har 4 leading digits accuracy to 11?

#define CATCH_CONFIG_RUNNER
#include "catch.hpp"
#include <iostream>
#include <cmath>
#include <array>
#include <vector>
#include <fstream>
#include <armadillo>
//#include "jacobi.h"

using namespace std;
using namespace arma;

//Function prototypes:
void printMatrixA(double** A, int N);
double maxoffdiag(double** A, int& k, int& l, int N);
vector<int> jacobi_method(double** A, double** R, int N);  //returnerer en vektor med ints
void rotate ( double ** A, double ** R, int k, int l, int N );

//Unit test (skal sjekke om jacobi finner rett egenverdier):
TEST_CASE("5x5 egenverdi-test") {
    int n = 5;
    double tolerance = pow(10.0, -2.0);
    double eigenvalue0_known = 1;
    double eigenvalue0;
    vector <int> idx(n);
    double eigval;
    vector <double> eigvec;

    //Set up 5x5 identity matrix:
    double** testMatrix_1 = new double*[n];
    for (int i = 0; i < n; i++){
        testMatrix_1[i] = new double[n];
    }
    for (int i = 0; i < n; i++){
        for (int j = 0; j < n; j++){
            if (i == j){
                testMatrix_1[i][j] = 1.0;
            }
            else{
                testMatrix_1[i][j] = 0.0;
            }
        }
    }

    double** R = new double*[n];
    for (int i = 0; i < n; i++){
        R[i] = new double[n];
    }

    // solve with jacobi:
    idx = jacobi_method(testMatrix_1, R, n);  //idx er vektor som inneholder indeksene til de 3 laveste egenverdiene
    eigenvalue0 = testMatrix_1[idx[0]][idx[0]];  //Den laveste egenverdien
    REQUIRE(abs(eigenvalue0 - eigenvalue0_known) < tolerance); //Krever at eigenvalue0 skal vaere naerme nok eigenvalue0_known

}//Slutt test

//Unit test (skal sjekke om maxoffdiag faktisk finner det storste offdiag elementet):
TEST_CASE("5x5 maks offdiag element test") {
    int n = 5;
    double tolerance = pow(10.0, -2.0);
    double maks_known = 0;
    int k;
    int l;

    //Set up 5x5 matrix with:
    double** testMatrix_2 = new double*[n];
    for (int i = 0; i < n; i++){
        testMatrix_2[i] = new double[n];
    }
    for (int i = 0; i < n; i++){
        for (int j = 0; j < n; j++){
            if (i == j){
                testMatrix_2[i][j] = 1.0;
            }
            else{
                testMatrix_2[i][j] = 0.0;
            }
        }

    }
    // use maxoffdiag:
    double maks = maxoffdiag(testMatrix_2, k, l, n);
    cout << "maks offdiag element: " << maks << endl;

    REQUIRE(abs(maks - maks_known) < tolerance);
}//Slutt test

int main(int argc, char* argv[]){

    int result = Catch::Session().run(argc, argv);  //kjorer testene

    //Definerer var som trengs for a lage matrisen A
    //N = 200 er passe?
    int N = 200;
    double rho_0 = 0;
    double rho_max = 8;
    double h = (rho_max - rho_0)/((double)N);
    double hh = h*h;
    double e = -1/hh;
    double w_r = 0.25;
    //double w_r = 5;
    double w_r2 = w_r*w_r;

    //Lager vektorer som trengs:
    double* rho = new double[N];
    double* V = new double[N];
    double* d = new double[N];
    double* eigenvec0 = new double[N];
    //double* eigenvalues = new double[N];
    vec eigval; // armadillo vektor

    //Fyller vektorene
    for (int i = 0; i < N; i++){
        rho[i] = rho_0 + (i+1)*h;
        V[i] = w_r2*rho[i]*rho[i] + 1.0/rho[i];
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

    //Maa definere B som en armadillo-matrise for a kunne bruke eig_sym
    mat B(N, N);
    for (int i = 0; i < N; i++){   //Fyller inn B slik at den blir lik som A
        for (int j = 0; j < N; j++){
            if (i == j-1){
                B(i,j) = e;
            } else if (i == j+1){
                B(i,j) = e;
            } else {
                B(i,j) = 0.0;
            }
        }
    }

    for (int i = 0; i < N; i++){
        B(i, i) = d[i];
    }

    eig_sym(eigval, B);//Bruker funksjonen eig_sym som lagrer egenverdiene i vektoren eigval
    eigval = sort(eigval); //Sorterer elementene i eigval fra minst til storst
    eigval(span(0,2)).print("De tre laveste egenverdiene fra eigsym:"); //Skriver ut elementene 0, 1 og 2 i eigval
    //eigval.print("Armadillo eigval: ");  //Skriver ut alle egenverdiene i eigval

    vector <int> ind = jacobi_method(A, R, N); //Kaller jacobi_method. ind inneholder naa indeksene til egenverdiene

    //Skriver til filen egenverdier.dat
    ofstream outFile;
    outFile.open("egenverdier.dat", std::ios::out);

    for (int j = 0; j < N; j++){
        eigenvec0[j] = R[ind[0]][j];  //eigenvec0 blir egenvektoren til grunntilstanden (dvs til laveste egenverdi)
        outFile << eigenvec0[j] << endl;  //Skriver egenvektoren til fil
    }
    outFile.close();

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
vector<int> jacobi_method(double** A, double** R, int N) {

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

    vector<double>diag(N);  //Deklarerer vektor diag med lengde N

    for (int i = 0; i < N; i++){
        diag[i] = A[i][i];  //Fyller diag med A sine diagonalelementer
    }

    //finner indeks til minste egenverdiene og lagrer indeksen i lowestIndex
    int ll = 3;  //ll (3) er antallet egenverdier jeg bryr meg om
    vector<double> lowest(ll);
    vector<int> lowestIndex(ll);  //vektor som skal inneholde indeksene til de minste egenverdiene
    for (int i = 0; i < ll; i++) {
        int index = 0;
        double low = 200;

        for (unsigned int j = 0; j < diag.size(); j++) {
            if (diag[j] < low) {
                low = diag[j];
                index = j;
            }
        }
        lowest[i] = low;
        lowestIndex[i] = index;
        diag[index] = 91385938659368;
    }

    //printMatrixA(A,N);

    cout << "De tre laveste egenverdiene fra Jacobi-metoden:" << endl;
    for (int i = 0; i < ll; i++) {
        cout << lowest[i] << " tilhorende index: " << lowestIndex[i] << endl;
    }

    sort(diag.begin(), diag.end());  //Sorterer vektoren fra minst til storste verdier
    //Skriver ut diag sine 3 forste elementer (3 forste egenverdiene):
    //cout << "3 forste egenverdier: " << diag[lowestIndex[0]] << " " << diag[1] << " " << diag[2] << endl;
    //cout << "3 forste egenverdier: " << diag[0] << " " << diag[1] << " " << diag[2] << endl;


    return lowestIndex;  //Inneholder indeksene til de 3 laveste egenverdiene
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
