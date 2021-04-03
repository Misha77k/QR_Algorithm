#include<stdlib.h>
#include<stdio.h>
#include<cstdlib>
#include<math.h>
#include<time.h>
#define EPS 10e-50//точность
double func(int i,int j);
int ReadfromFile(char *fp,double *A,int N);
double vect_norma(double *a,int n);
double scalar_multi(double *a,double *b,int n);
double matr_norma(double *A,int N);
int Hausholder_matr(double *A,int N,int i,double *c,double *w);
int Givens_matr(double *A,int N,int i,int M,double *s,double *r,double *p);
void QR(double *A,double *x,int N,double *s,double *r,double *p);
double nevyazka_tr(double *A,double *x,int N);
double nevyazka_len(double *A,double *x,int N);
void print(double *A,int N);
void print_part(double *A,int N,int K);
