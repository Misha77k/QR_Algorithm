#include<stdlib.h>
#include<stdio.h>
#include<cstdlib>
#include<math.h>
#include<time.h>
#include"QR.h"
#define EPS 10e-50//точность
#define FPS 10e-10
double func(int i,int j){
    return abs(i-j);
}
int ReadfromFile(char*fp,double *A,int N){
    FILE *f;
    int K=0;
    // int M=0;
    char a;
    f=fopen(fp,"r");
    if(f==NULL){
        printf("файл не открылся\n");
        return -1;
    }
    for(int i=0;i<N;i++)
        for(int j=0;j<N;j++){
            if(fscanf(f,"%lf",&A[i*N+j])!=1){
                printf("error\n");
                return -1;
            }
            K++;
        }
    int q=fscanf(f,"%c",&a);
    if (fscanf(f,"%c",&a)==1){
        printf("Несовпадение размерностей\n");
        fclose(f);
        return -1;
    }
    if(K!=N*N){
        printf("Несовпадение размерностей\n");
        fclose(f);
        return -1;
    }
    fclose(f);
    return 0;
}
double vect_norma(double *a,int n){
    double d=0.0;
    for(int i=0;i<n;i++)
        d=d+a[i]*a[i];
    return d;
}
double scalar_multi(double *a,double *b,int n){
    double d=0.0;
    for(int i=0;i<n;i++)
        d=d+a[i]*b[i];
    return d;
}
double matr_norma(double *A,int N){
    double max=0.0;
    double cur;
    for(int i=0;i<N;i++){
        cur=0.0;
        for(int j=0;j<N;j++)
            cur=cur+fabs(A[i*N+j]);
        if(cur>max)
            max=cur;
    }
    return max;
}
int Hausholder_matr(double *A,int N,int i,double *c,double *w){
    for(int j=0;j<N-i-1;j++)
        c[j]=A[(j+i+1)*N+i];
    double q,g;
    g=sqrt(vect_norma(c,N-i-1));
    c[0]=c[0]+g;
    q=vect_norma(c,N-i-1);
    // if(fabs(q)<FPS)
    //   return -1;
    double y;
    y=2/q;
    for(int j=0;j<N-i-1;j++){
        w[j]=0.0;
        for(int k=0;k<N-i-1;k++)
            w[j]=w[j]+A[(j+i+1)*N+k+i+1]*c[k];
    }
    
    for(int j=0;j<N-i-1;j++)
        w[j]=w[j]*y;
    double x=0.0;
    x=scalar_multi(c,w,N-i-1);
    x=x*y;
    for(int j=0;j<N-i-1;j++){
        for(int k=0;k<N-i-1;k++){
            A[(j+i+1)*N+k+i+1]=A[(j+i+1)*N+k+i+1]-c[j]*w[k]-w[j]*c[k]+x*c[j]*c[k];
        }
    }
    A[i*N+i+1]=g;
    A[(i+1)*N+i]=g;
    for(int j=i+2;j<N;j++){
        A[j*N+i]=0;
        A[i*N+j]=0;
    }
    return 0;
}
int Givens_matr(double *A,int N,int i,int M,double *s,double *r,double *p){
    //for(int j=0;j<2;j++)
    //  s[j]=A[(i+j)*M+i];
    s[0]=A[i*M+i];
    s[1]=A[(i+1)*M+i];
    double q=0.0;
    q=sqrt(vect_norma(s,2));
    if(q==0)
        return -1;
    s[0]=s[0]/q;
    s[1]=-s[1]/q;
    /*for(int j=0;j<N;j++){
     r[j]=0;
     p[j]=0;
     }*/
    /*for(int j=0;j<N;j++){
     if(j==i){
     r[j]=s[0];
     r[j+1]=-s[1];
     p[j]=s[1];
     p[j+1]=s[0];
     }
     }*/
    r[i]=s[0];
    //r[i+1]=-s[1];
    p[i]=s[1];
    //p[i+1]=s[0];
    //r[i+1]=-p[i];
    //p[i+1]=r[i];
    double y;
    double z;
    for(int j=0;j<N;j++){
        y=0.0;
        z=0.0;
        /* for(int k=0;k<N;k++){
         y=y+r[k]*A[k*M+j];
         z=z+p[k]*A[k*M+j];
         }*/
        y=y+r[i]*A[i*M+j]-p[i]*A[(i+1)*M+j];
        z=z+p[i]*A[i*M+j]+r[i]*A[(i+1)*M+j];
        A[i*M+j]=y;
        A[(i+1)*M+j]=z;
        //y=0.0;
        //z=0.0;
        /*for(int k=0;k<N;k++){
         y=y+r[k]*Q[j*M+k];
         z=z+p[k]*Q[j*M+k];
         }*/
        /*y=y+r[i]*Q[j*M+i]+r[i+1]*Q[j*M+i+1];
         z=z+p[i]*Q[j*M+i]+p[i+1]*Q[j*M+i+1];
         Q[j*M+i]=y;
         Q[j*M+i+1]=z;*/
    }
    return 0;
}
void QR(double *A,double *x,int N,double *s,double *r,double *p){
    //int k=0;
    //clock_t t1;
    int l=0;
    double nor=matr_norma(A,N);
    double q;
    //int u=0;
    int M=N;
    while(N!=2){
        q=A[(N-1)*M+N-1];
        for(int i=0;i<N;i++)
            A[i*M+i]=A[i*M+i]-q;
        /*for(int i=0;i<N;i++)
         for(int j=0;j<N;j++){
         if(i==j)
         Q[i*M+j]=1;
         else Q[i*M+j]=0;
         }*/
        //t1=clock();
        for(int i=0;i<N-1;i++){
            l=Givens_matr(A,N,i,M,s,r,p);
        }
        /*printf("\n");
         for(int i=0;i<N;i++){
         for(int j=0;j<N;j++)
         printf("%lf ",A[i*N+j]);
         printf("\n");
         }
         printf("\n");
         for(int i=0;i<N-1;i++)
         printf("rrr=%lf ",r[i]);
         printf("\n");
         for(int i=0;i<N-1;i++)
         printf("ppp=%lf ",p[i]);
         printf("\n");*/
        // t1=clock()-t1;
        // printf("  time1=  %lf\n",((double)t1)/CLOCKS_PER_SEC);
        /*for(int i=0;i<N;i++)
         for(int j=i;j<N;j++)
         R[i*M+j]=A[i*M+j];*/
        /* for(int i=0;i<N;i++)
         for(int j=0;j<N;j++)
         A[i*M+j]=0.0;*/
        //multi_matr(R,Q,A,N,M);
        double y;
        double z;
        for(int i=0;i<N;i++)
            for(int j=0;j<N-1;j++){
                y=0.0;
                z=0.0;
                y=y+A[i*M+j]*r[j]-A[i*M+j+1]*p[j];
                z=z+A[i*M+j]*p[j]+A[i*M+j+1]*r[j];
                A[i*M+j]=y;
                A[i*M+j+1]=z;
            }
        /* for(int i=0;i<N;i++){
         for(int j=0;j<N;j++)
         printf("%lf ",A[i*N+j]);
         printf("\n");
         }
         printf("\n");*/
        for(int i=0;i<N;i++)
            A[i*M+i]=A[i*M+i]+q;
        if(nor==0){
            x[N-1]=0;
            N--;
        }
        if(fabs(A[(N-1)*M+N-2])<EPS*nor){
            x[N-1]=q;
            N--;
            // u++;
        }
    }
    double D;
    D=(A[0]+A[M+1])*(A[0]+A[M+1])-4*(A[0]*A[M+1]-A[M]*A[1]);
    x[0]=(A[M+1]+A[0]+sqrt(D))/2;
    x[1]=(A[M+1]+A[0]-sqrt(D))/2;
}
double nevyazka_tr(double *A,double *x,int N){
    double trA=0.0;
    double trx=0.0;
    for(int i=0;i<N;i++){
        trA=trA+A[i*N+i];
        trx=trx+x[i];
    }
    return(fabs(trA-trx));
}
double nevyazka_len(double *A,double *x,int N){
    double lenA=0.0;
    double lenx=0.0;
    for(int i=0;i<N;i++)
        for(int j=0;j<N;j++)
            lenA=lenA+A[i*N+j]*A[i*N+j];
    for(int i=0;i<N;i++)
        lenx=lenx+x[i]*x[i];
    return(fabs(lenA-lenx));
}
void print(double *A,int N){
    for(int i=0;i<N;i++){
        for(int j=0;j<N;j++)
            printf("%lf  ",A[i*N+j]);
        printf("\n");
    }
}
void print_part(double *A,int N,int K){
    if(K<N){
        for(int i=0;i<K;i++){
            for(int j=0;j<K;j++)
                printf("%lf ",A[i*N+j]);
            printf("\n");
        }
    }
    if(K>=N)
        print(A,N);
}


