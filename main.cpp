#include <iostream>
#include <cstdlib>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <string>
#include <armadillo>
#include "time.h"

//2.d Quantum dots in three dimensions

using namespace std;
using namespace arma;

//The functions below the main function
arma::mat make_tridiag_A(int n);
void maxval(mat A, double max, int &k, int &l, int n);
void Jacobi(mat &A, mat R, int k, int l, int n);

int main (int argc, char* argv[])
{
  //Choosing n from the command line
  int n = atoi(argv[1]);

  //The eigenvectors for the starting diagonal matrice is an identity matrice
  arma::mat R = eye(n,n);
  arma::mat A = make_tridiag_A(n);

  //The eigenvals for A using armadillo
  arma::vec eigvals = arma::eig_sym(A);
  cout << "The analytical eigenvalues are " << endl;
  //Analytival Eigenvalues
  //double analytic_eig[n];
  vec analytic_eig = {3,5,7,11,13,17,19,23,29,31};
  for(int i=0 ; i < n ; i++) {
      cout << analytic_eig[i] << endl;
  }

  double tolerance = 1.0E-10;
  double max=3;
  int k, l;
  int count = 0;
  clock_t start, finish;
  start = clock();
  while (abs(max) > tolerance, count < atoi(argv[2])){
     maxval(A, max, k, l, n);
     Jacobi(A, R, k, l, n);
     count = count + 1;
  }
  finish = clock();
  double timeused = (double) (finish - start)/(CLOCKS_PER_SEC );

  cout << "The armadillo eigenvalues are" << endl;
  cout << eigvals << endl;
  cout << setiosflags(ios::showpoint | ios::uppercase);
  cout << "n = " << n << endl;
  cout << setw(20) << "The time used = " << timeused  << endl;
  cout << setw(20) << "Number of iterations in the while loop = "<< count << endl;
  cout << "The tolerance factor is " << tolerance << endl;

  cout << "The new matrice A: " << endl;
  cout <<  A << endl;
  cout << "The eigenvector matrice R:" << endl;
  cout << R << endl;

}


arma::mat make_tridiag_A(int n){
    //Making the first A matrice, now with added potential on the diagonal elements
    double h = 1/(double) (n+1); //stepsize
    double d = 2.0/(h*h);
    double a = -1.0/(h*h);


    arma::mat A = arma::mat(n, n);

    for(int i=0 ; i < n ; i++) {
       for(int j=0 ; j < n ; j++) {
           double h = 1/(double) (n+1);
           double p[n];
           double V[n];
           p[i] = i*h;
           V[i] = p[i]*p[i]; //Harmonic Oscillator potential

           if (i == j){
               A(i,j) = d + V[i];
           }
           else if (abs(i-j)<=1) {
               A(i,j) = a;

           }
           else {
               A(i,j) = 0.0;
           }

       }
     }

    return A;
}



void maxval(mat A, double max, int &k, int &l, int n){
    //Finding the maximal off-diagonal value
    for (int i = 0; i < n; ++i)
    {
        for (int j = i+1; j < n; ++j)
        {
            double maxoffdiag;
            double aij = fabs(A(i,j));
            if (aij > maxoffdiag){
                maxoffdiag = aij;
                k = i;
                l = j;
            }

         }
    }

}

void Jacobi(mat &A, mat R, int k, int l, int n){
    //Jacobi transformation
    double tau = (A(l,l)-A(k,k))/(2*A(k,l));
    double t = -tau + sqrt(1+tau*tau);
    double c = 1/sqrt(1+t*t);
    double s = t*c;
    double akk = A(k,k);
    double all = A(l,l);
    A(k,l) = 0.0;
    A(l,k) = 0.0;
    A(k,k) = akk*c*c-2.0*c*s*A(k,l)+all*s*s;
    A(l,l) = all*c*c+2.0*c*s*A(k,l)+akk*s*s;
    for (int i = 0; i < n; i++){
        if (i !=k && i !=l){
            double aik = A(i,k);
            double ail = A(i,l);
            A(i,k) = aik*c-ail*s;
            A(k,i) = A(i,k);
            A(i,l) = ail*c+aik*s;
            A(l,i) = A(i,l);
            //Making the new eigenvalu
            double rik = R(i,k);
            double ril = R(i,l);
            R(i,k) = c*rik-s*ril;
            R(i,l) = c*ril+s*rik;

        }
    }

}



