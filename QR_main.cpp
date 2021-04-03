#include<stdlib.h>
#include<stdio.h>
#include<cstdlib>
#include<math.h>
#include<time.h>
#include"QR.h"
#define EPS 10e-50//точность
int main(int argc,char **argv){
    double *A;
    double *B;
    // double *E;
    //double *Q;
    //double *R;
    double *x;
    int N;
    clock_t t;
    if(argc!=2 && argc!=3){
        printf("error\n");
        return -1;
    }
    else{
        N=atoi(argv[1]);
    }
    if(N<=0){
        printf("Неверная размерность\n");
        return -1;
    }
    // printf("Введи размерность матрицы\n");
    //scanf("%d",&N);
    A=(double *)malloc(N*N*sizeof(double));
    //E=(double *)malloc(N*N*sizeof(double));
    //Q=(double *)malloc(N*N*sizeof(double));
    //R=(double *)malloc(N*N*sizeof(double));
    
    /*for(int i=0;i<N;i++)
     for(int j=0;j<N;j++){
     if(i==j)
     Q[i*N+j]=1;
     else Q[i*N+j]=0;
     }*/
    //for(int i=0;i<N*N;i++)
    //  Q[i]=E[i];
    /*for(int i=0;i<N*N;i++)
     if(scanf("%lf",&A[i])!=1)
     return -1;*/
    if(argc==2){
        for(int i=0;i<N;i++)
            for(int j=0;j<N;j++)
                A[i*N+j]=func(i,j);
    }
    if(argc==3){
        int k = ReadfromFile(argv[2],A,N);
        if(k==-1){
            free(A);
            //free(Q);
            //free(R);
            return -1;
        }
    }
    /*for(int i=0;i<N;i++)
     for(int j=0;j<i;j++)
     R[i*N+j]=0;*/
    /*double trA=0.0;
     for(int i=0;i<N;i++)
     trA=trA+A[i*N+i];*/
    B=(double *)malloc(N*N*sizeof(double));
    for(int i=0;i<N;i++)
        for(int j=0;j<N;j++)
            B[i*N+j]=A[i*N+j];
    x=(double *)malloc(N*sizeof(double));
    if(N==1){
        x[0]=A[0];
        printf("%lf",x[0]);
        free(x);
        free(A);
        free(B);
        return 0;
    }
    else{
        for(int i=0;i<N;i++)
            x[i]=0.0;
        double *c;
        double *w;
        double *s;
        double *p;
        double *r;
        p=(double *)malloc(N*sizeof(double));
        r=(double *)malloc(N*sizeof(double));
        s=(double *)malloc(2*sizeof(double));
        int k=0;
        //t=clock();
        for(int i=0;i<N-2;i++){
            c=(double *)malloc((N-i-1)*sizeof(double));
            w=(double *)malloc((N-i-1)*sizeof(double));
            k=Hausholder_matr(A,N,i,c,w);
            free(c);
            free(w);
        }
        if(k==-1)
            return -1;
        t=clock();
        QR(A,x,N,s,r,p);
        t=clock()-t;
        printf("  time=  %lf\n",((double)t)/CLOCKS_PER_SEC);
        /*double trx=0.0;
         for(int i=0;i<N;i++)
         trx=trx+x[i];*/
        /*for(int i=0;i<N;i++){
         for(int j=0;j<N;j++)
         printf("%lf ",A[i*N+j]);
         printf("\n");
         }*/
        printf("\n");
        printf("{");
        for(int i=0;i<N;i++)
            printf("%lf  ",x[i]);
        printf("}");
        printf("\n");
        double nev1;
        nev1=nevyazka_tr(B,x,N);
        printf("невязка_след=%le\n",nev1);
        double nev2;
        nev2=nevyazka_len(B,x,N);
        printf("невязка_длина=%le\n",nev2);
        free(p);
        free(r);
        free(s);
        free(x);
        //free(R);
        //free(E);
        //free(Q);
        free(B);
        free(A);
        return 0;
    }
}
