#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "lapack_prototypes.h"
#include "Generate_MLWF.h"

void lapack_dstevx2(INTEGER N, INTEGER EVmax, double *D, double *E, double *W, dcomplex **ev, int ev_flag)
{
    int i,j;

    char  *JOBZN="N";
    char  *JOBZV="V";
    char  *RANGE="I";

    double VL,VU; /* dummy */
    INTEGER IL,IU; 
    double ABSTOL=LAPACK_ABSTOL;
    INTEGER M;
    double *Z;
    INTEGER LDZ;
    double *WORK;
    INTEGER *IWORK;
    INTEGER *IFAIL;
    INTEGER INFO;

    IL = 1;
    IU = EVmax;

    M = IU - IL + 1;
    LDZ = N;

    Z = (double*)malloc(sizeof(double)*LDZ*N);
    WORK = (double*)malloc(sizeof(double)*5*N);
    IWORK = (INTEGER*)malloc(sizeof(INTEGER)*5*N);
    IFAIL = (INTEGER*)malloc(sizeof(INTEGER)*N);

    if (ev_flag){
        F77_NAME(dstevx,DSTEVX)( JOBZV, RANGE, &N,  D, E, &VL, &VU, &IL, &IU, &ABSTOL,
                &M, W, Z, &LDZ, WORK, IWORK, IFAIL, &INFO );
    }
    else{
        F77_NAME(dstevx,DSTEVX)( JOBZN, RANGE, &N,  D, E, &VL, &VU, &IL, &IU, &ABSTOL,
                &M, W, Z, &LDZ, WORK, IWORK, IFAIL, &INFO );
    }

    /* store eigenvectors */

    if (ev_flag==1){
        for (i=0; i<EVmax; i++) {
            for (j=0; j<N; j++) {
                ev[i+1][j+1].r = Z[i*N+j];
                ev[i+1][j+1].i = 0.0;
            }
        }
    }

    for (i=EVmax; i>=1; i--){
        W[i]= W[i-1];
    }

    if (INFO>0) {
    }

    if (INFO<0) {
        printf("info=%d in dstevx_\n",INFO);
        exit(0);
    }

    free(Z);
    free(WORK);
    free(IWORK);
    free(IFAIL);
}



int Eigen_zheevx(dcomplex **A, double *W, int N0, int MaxN, int ev_flag)
{
    char *JOBZ0="N";
    char *JOBZ1="V";
    char *UPLO="L";
    char *RANGE="I";
    INTEGER N,LDA,LDZ;
    double VL,VU,sum;
    INTEGER  IL,IU;
    double ABSTOL=LAPACK_ABSTOL;
    INTEGER M;
    dcomplex *Z;
    dcomplex *WORK;
    INTEGER LWORK;
    double *RWORK;
    INTEGER  *IWORK;
    INTEGER *IFAIL,INFO;
    dcomplex *A0;
    int i,j,k,po;

    N = N0;
    LDA = N;
    LDZ = N;

    IL = 1;
    IU = MaxN; 

    A0=(dcomplex*)malloc(sizeof(dcomplex)*N*N);
    for (i=0; i<N*N; i++){
        A0[i].r = 0.0;
        A0[i].i = 0.0;
    }

    for (i=1;i<=N;i++) {
        for (j=1;j<=N;j++) {
            A0[(j-1)*N+i-1].r = A[i][j].r;
            A0[(j-1)*N+i-1].i = A[i][j].i;
        }
    }

    Z=(dcomplex*)malloc(sizeof(dcomplex)*N*N);
    for (i=0; i<N*N; i++){
        Z[i].r = 0.0;
        Z[i].i = 0.0;
    }

    LWORK= 3*N;
    WORK=(dcomplex*)malloc(sizeof(dcomplex)*LWORK);
    for (i=0; i<LWORK; i++){
        WORK[i].r =0.0;
        WORK[i].i =0.0;
    }

    RWORK=(double*)malloc(sizeof(double)*7*N);
    for (i=0; i<7*N; i++) RWORK[i] = 0.0;

    IWORK=(INTEGER*)malloc(sizeof(INTEGER)*5*N);
    for (i=0; i<5*N; i++) IWORK[i] = 0;

    IFAIL=(INTEGER*)malloc(sizeof(INTEGER)*N);
    for (i=0; i<N; i++) IFAIL[i] = 0;


    F77_NAME(zheevx,ZHEEVX)(JOBZ1, RANGE, UPLO, &N, A0, &LDA, &VL, &VU, &IL, &IU,
                            &ABSTOL, &M, W, Z, &LDZ, WORK, &LWORK, RWORK,
                            IWORK, IFAIL, &INFO );

    if (ev_flag==1){

        for (i=1; i<=N; i++) {
            for (j=1; j<=N; j++) {
                A[i][j].r = Z[(j-1)*N+i-1].r;
                A[i][j].i = Z[(j-1)*N+i-1].i;
            }
        }
    }

    for (i=N; i>=1; i--) {
        W[i] = W[i-1];
    }

    free(IFAIL); free(IWORK); free(RWORK); free(WORK); free(Z); free(A0); 

    return INFO;
}


#pragma optimization_level 1
void Eigen_HH(dcomplex **ac, double *ko, int n, int EVmax, int ev_flag)
{

    double ABSTOL=LAPACK_ABSTOL;

    dcomplex **ad,*u,*b1,*p,*q,tmp0,tmp1,tmp2,tmp3,u1,u2,p1;
    dcomplex ss0,ss1,ss2,ss3,ss,p10,p11,p12,p13;
    double *D,*E,*uu,*alphar,*alphai,
            s1,s2,s3,r,
            sum,ar,ai,br,bi,e,
            a1,a2,a3,a4,a5,a6,b7,r0,r1,r2,
            r3,x1,x2,xap,
            bb,bb1,ui,uj,uij;

    int jj,jj1,jj2,k,ii,ll,i3,i2,j2,
        i,j,i1,j1,n1,n2,ik,ks,i1s,
        jk,po1,nn,count;

    double Stime, Etime;
    double Stime1, Etime1;
    double Stime2, Etime2;
    double time1,time2;


    n2 = n + 5;

    ad = (dcomplex**)malloc(sizeof(dcomplex*)*n2);
    for (i=0; i<n2; i++){
        ad[i] = (dcomplex*)malloc(sizeof(dcomplex)*n2);
    }

    b1 = (dcomplex*)malloc(sizeof(dcomplex)*n2);
    u = (dcomplex*)malloc(sizeof(dcomplex)*n2);
    uu = (double*)malloc(sizeof(double)*n2);
    p = (dcomplex*)malloc(sizeof(dcomplex)*n2);
    q = (dcomplex*)malloc(sizeof(dcomplex)*n2);

    D = (double*)malloc(sizeof(double)*n2);
    E = (double*)malloc(sizeof(double)*n2);

    alphar = (double*)malloc(sizeof(double)*n2);
    alphai = (double*)malloc(sizeof(double)*n2);

    for (i=1; i<=(n+2); i++){
        uu[i] = 0.0;
    }


    for (i=1; i<=(n-1); i++){

        s1 = ac[i+1][i].r * ac[i+1][i].r + ac[i+1][i].i * ac[i+1][i].i;
        s2 = 0.0;

        u[i+1].r = ac[i+1][i].r;
        u[i+1].i = ac[i+1][i].i;
    
        for (i1=i+2; i1<=(n-1); i1+=2){
            s2 += ac[i1+0][i].r*ac[i1+0][i].r + ac[i1+0][i].i*ac[i1+0][i].i + ac[i1+1][i].r*ac[i1+1][i].r + ac[i1+1][i].i*ac[i1+1][i].i;
            u[i1+0] = ac[i1+0][i];
            u[i1+1] = ac[i1+1][i];
        }

        i1s = n + 1 - (n+1-(i+2))%2;

        for (i1=i1s; ((i+2)<=i1 && i1<=n); i1++){
            s2 += ac[i1][i].r*ac[i1][i].r + ac[i1][i].i*ac[i1][i].i;
            u[i1] = ac[i1][i];
        }
        s3 = fabs(s1 + s2);

        if ( ABSTOL<(fabs(ac[i+1][i].r)+fabs(ac[i+1][i].i)) ){
            if (ac[i+1][i].r<0.0)  s3 =  sqrt(s3);
            else                   s3 = -sqrt(s3);
        }else{
            s3 = sqrt(s3);
        }

        if ( ABSTOL<fabs(s2) || 1.0e-10<fabs(u[i+1].i) || i==(n-1) ){

            ss.r = ac[i+1][i].r;
            ss.i = ac[i+1][i].i;

            ac[i+1][i].r = s3;
            ac[i+1][i].i = 0.0;
            ac[i][i+1].r = s3;
            ac[i][i+1].i = 0.0;

            u[i+1].r = u[i+1].r - s3;
            u[i+1].i = u[i+1].i;
            
            u1.r = s3 * s3 - ss.r * s3;
            u1.i =         - ss.i * s3;
            u2.r = 2.0 * u1.r;
            u2.i = 2.0 * u1.i;
            
            e = u2.r/(u1.r*u1.r + u1.i*u1.i);
            ar = e*u1.r;
            ai = e*u1.i;

            alphar[i] = ar;
            alphai[i] = ai;

            uu[i] = u2.r;

            b1[i].r = ss.r - s3;
            b1[i].i = ss.i;

            r = 0.0;
            for (i1=i+1; i1<=n; i1++){

                p1.r = 0.0;
                p1.i = 0.0;
                for (j=i+1; j<=n; j++){
                    p1.r += ac[i1][j].r * u[j].r - ac[i1][j].i * u[j].i;
                    p1.i += ac[i1][j].r * u[j].i + ac[i1][j].i * u[j].r;
                }
                p[i1].r = p1.r / u1.r;
                p[i1].i = p1.i / u1.r;

                r += u[i1].r * p[i1].r + u[i1].i * p[i1].i;
            }
            r = 0.5*r / u2.r;

            br =  ar*r;
            bi = -ai*r;

            for (i1=i+1; i1<=n; i1++){
                tmp1.r = 0.5*(p[i1].r - (br * u[i1].r - bi*u[i1].i));
                tmp1.i = 0.5*(p[i1].i - (br * u[i1].i + bi*u[i1].r));
                q[i1].r = ar * tmp1.r - ai * tmp1.i; 
                q[i1].i = ar * tmp1.i + ai * tmp1.r; 
            }

            for (i1=i+1; i1<=n; i1++){
                tmp1.r = u[i1].r;
                tmp1.i = u[i1].i;
                tmp2.r = q[i1].r; 
                tmp2.i = q[i1].i; 
                for (j1=i+1; j1<=n; j1++){
                    ac[i1][j1].r -= ( tmp1.r * q[j1].r + tmp1.i * q[j1].i+tmp2.r * u[j1].r + tmp2.i * u[j1].i );
                    ac[i1][j1].i -= (-tmp1.r * q[j1].i + tmp1.i * q[j1].r-tmp2.r * u[j1].i + tmp2.i * u[j1].r );
                }
            }
        }
    }

    for (i=1; i<=n; i++){
        for (j=1; j<=n; j++){
            ad[i][j].r = ac[i][j].r;
            ad[i][j].i = ac[i][j].i;
        }
    }



    for (i=1; i<=n; i++){
        D[i-1] = ad[i][i  ].r;
        E[i-1] = ad[i][i+1].r;
    }

    /*
    if      (dste_flag==0) lapack_dstegr2(n,EVmax,D,E,ko,ac);
    else if (dste_flag==1) lapack_dstedc2(n,D,E,ko,ac);
    else if (dste_flag==2) lapack_dstevx2(n,EVmax,D,E,ko,ac,ev_flag);
    */


    lapack_dstevx2(n,EVmax,D,E,ko,ac,ev_flag);


    if (ev_flag==1){

        for (i=2; i<=n; i++){
            ad[i-1][i].r = b1[i-1].r;
            ad[i-1][i].i =-b1[i-1].i;
            ad[i][i-1].r = b1[i-1].r;
            ad[i][i-1].i = b1[i-1].i;
        }


        for (k=1; k<=(EVmax-3); k+=4){
            for (nn=1; nn<=n-1; nn++){
                if ( (1.0e-3*ABSTOL)<fabs(uu[n-nn])){

                    tmp0.r = 0.0;	tmp0.i = 0.0;
                    tmp1.r = 0.0;	tmp1.i = 0.0;
                    tmp2.r = 0.0;	tmp2.i = 0.0;
                    tmp3.r = 0.0;	tmp3.i = 0.0;

                    for (i=n-nn+1; i<=n; i++){
                        tmp0.r += ad[n-nn][i].r * ac[k+0][i].r - ad[n-nn][i].i * ac[k+0][i].i;
                        tmp0.i += ad[n-nn][i].i * ac[k+0][i].r + ad[n-nn][i].r * ac[k+0][i].i;

                        tmp1.r += ad[n-nn][i].r * ac[k+1][i].r - ad[n-nn][i].i * ac[k+1][i].i;
                        tmp1.i += ad[n-nn][i].i * ac[k+1][i].r + ad[n-nn][i].r * ac[k+1][i].i;

                        tmp2.r += ad[n-nn][i].r * ac[k+2][i].r - ad[n-nn][i].i * ac[k+2][i].i;
                        tmp2.i += ad[n-nn][i].i * ac[k+2][i].r + ad[n-nn][i].r * ac[k+2][i].i;

                        tmp3.r += ad[n-nn][i].r * ac[k+3][i].r - ad[n-nn][i].i * ac[k+3][i].i;
                        tmp3.i += ad[n-nn][i].i * ac[k+3][i].r + ad[n-nn][i].r * ac[k+3][i].i;
                    }

                    ss0.r = (alphar[n-nn]*tmp0.r - alphai[n-nn]*tmp0.i) / uu[n-nn];
                    ss0.i = (alphar[n-nn]*tmp0.i + alphai[n-nn]*tmp0.r) / uu[n-nn];

                    ss1.r = (alphar[n-nn]*tmp1.r - alphai[n-nn]*tmp1.i) / uu[n-nn];
                    ss1.i = (alphar[n-nn]*tmp1.i + alphai[n-nn]*tmp1.r) / uu[n-nn];

                    ss2.r = (alphar[n-nn]*tmp2.r - alphai[n-nn]*tmp2.i) / uu[n-nn];
                    ss2.i = (alphar[n-nn]*tmp2.i + alphai[n-nn]*tmp2.r) / uu[n-nn];

                    ss3.r = (alphar[n-nn]*tmp3.r - alphai[n-nn]*tmp3.i) / uu[n-nn];
                    ss3.i = (alphar[n-nn]*tmp3.i + alphai[n-nn]*tmp3.r) / uu[n-nn];

                    for (i=n-nn+1; i<=n; i++){
                        ac[k+0][i].r -= ss0.r * ad[n-nn][i].r + ss0.i * ad[n-nn][i].i;
                        ac[k+0][i].i -=-ss0.r * ad[n-nn][i].i + ss0.i * ad[n-nn][i].r;

                        ac[k+1][i].r -= ss1.r * ad[n-nn][i].r + ss1.i * ad[n-nn][i].i;
                        ac[k+1][i].i -=-ss1.r * ad[n-nn][i].i + ss1.i * ad[n-nn][i].r;

                        ac[k+2][i].r -= ss2.r * ad[n-nn][i].r + ss2.i * ad[n-nn][i].i;
                        ac[k+2][i].i -=-ss2.r * ad[n-nn][i].i + ss2.i * ad[n-nn][i].r;

                        ac[k+3][i].r -= ss3.r * ad[n-nn][i].r + ss3.i * ad[n-nn][i].i;
                        ac[k+3][i].i -=-ss3.r * ad[n-nn][i].i + ss3.i * ad[n-nn][i].r;
                    }
                }
            }
        }

        ks = EVmax - EVmax%4 + 1;

        for (k=ks; k<=EVmax; k++){
            for (nn=1; nn<=n-1; nn++){
                if ( (1.0e-3*ABSTOL)<fabs(uu[n-nn])){

                    tmp1.r = 0.0;
                    tmp1.i = 0.0;

                    for (i=n-nn+1; i<=n; i++){
                        tmp1.r += ad[n-nn][i].r * ac[k][i].r - ad[n-nn][i].i * ac[k][i].i;
                        tmp1.i += ad[n-nn][i].i * ac[k][i].r + ad[n-nn][i].r * ac[k][i].i;
                    }

                    ss.r = (alphar[n-nn]*tmp1.r - alphai[n-nn]*tmp1.i) / uu[n-nn];
                    ss.i = (alphar[n-nn]*tmp1.i + alphai[n-nn]*tmp1.r) / uu[n-nn];

                    for (i=n-nn+1; i<=n; i++){
                        ac[k][i].r -= ss.r * ad[n-nn][i].r + ss.i * ad[n-nn][i].i;
                        ac[k][i].i -=-ss.r * ad[n-nn][i].i + ss.i * ad[n-nn][i].r;
                    }
                }
            }
        }

    
        for (j=1; j<=EVmax; j++){
            sum = 0.0;
            for (i=1; i<=n; i++){
                sum += ac[j][i].r * ac[j][i].r + ac[j][i].i * ac[j][i].i;
            }
            sum = 1.0/sqrt(sum);
            for (i=1; i<=n; i++){
                ac[j][i].r = ac[j][i].r * sum;
                ac[j][i].i = ac[j][i].i * sum;
            }
        }

        for (i=1; i<=n; i++){
            for (j=(i+1); j<=n; j++){
                tmp1 = ac[i][j];
                tmp2 = ac[j][i];
                ac[i][j] = tmp2;
                ac[j][i] = tmp1;
            }
        }
    }


    for (i=0; i<n2; i++){
        free(ad[i]);
    }
    free(ad);

    free(b1);
    free(u);
    free(uu);
    free(p);
    free(q);
    free(D);
    free(E);
    free(alphar);
    free(alphai);
}


void EigenBand_lapack(dcomplex **A, double *W, int N0, int MaxN, int ev_flag)
{
    int info;

    info = Eigen_zheevx(A,W,N0,MaxN,ev_flag);
    if (info!=0){
        Eigen_HH(A,W,N0,MaxN,ev_flag);
    }
}



#pragma optimization_level 1
void Overlap_Band_Wannier(double ****OLP, dcomplex **S, int Natom, int *FNAN, int **natn, int **ncn, int **atv_ijk, int *Total_NumOrbs, int *MP, double k1, double k2, double k3)
{
    int i,j,wanA,wanB,tnoA,tnoB,Anum,Bnum,NUM,GA_AN,LB_AN,GB_AN;
    int l1,l2,l3,Rn,n2;
    double **S1,**S2;
    double kRn,si,co,s;

    Anum = 1;
    for (i=0; i<Natom; i++){
        MP[i] = Anum;
        Anum += Total_NumOrbs[i];
    }
    NUM = Anum - 1;

    /****************************************************
                         Allocation
    ****************************************************/

    n2 = NUM + 2;

    S1 = (double**)malloc(sizeof(double*)*n2);
    for (i=0; i<n2; i++){
        S1[i] = (double*)malloc(sizeof(double)*n2);
    }

    S2 = (double**)malloc(sizeof(double*)*n2);
    for (i=0; i<n2; i++){
        S2[i] = (double*)malloc(sizeof(double)*n2);
    }

    /****************************************************
                         set overlap
    ****************************************************/

    S[0][0].r = NUM;

    for (i=1; i<=NUM; i++){
        for (j=1; j<=NUM; j++){
            S1[i][j] = 0.0;
            S2[i][j] = 0.0;
        }
    }

    for (GA_AN=0; GA_AN<Natom; GA_AN++){
        tnoA = Total_NumOrbs[GA_AN];
        Anum = MP[GA_AN];

        for (LB_AN=0; LB_AN<=FNAN[GA_AN]; LB_AN++){
            GB_AN = natn[GA_AN][LB_AN]-1;
            Rn = ncn[GA_AN][LB_AN];
            tnoB = Total_NumOrbs[GB_AN];

            l1 = atv_ijk[Rn][0];
            l2 = atv_ijk[Rn][1];
            l3 = atv_ijk[Rn][2];
            kRn = k1*(double)l1 + k2*(double)l2 + k3*(double)l3;

            si = sin(2.0*PI*kRn);
            co = cos(2.0*PI*kRn);
            Bnum = MP[GB_AN];
            for (i=0; i<tnoA; i++){
                for (j=0; j<tnoB; j++){
                    s = OLP[GA_AN][LB_AN][i][j];
                    S1[Anum+i][Bnum+j] += s*co;
                    S2[Anum+i][Bnum+j] += s*si;
                }
            }
        }
    }

    for (i=1; i<=NUM; i++){
        for (j=1; j<=NUM; j++){
            S[i][j].r =  S1[i][j];
            S[i][j].i =  S2[i][j];
        }
    }

    /****************************************************
                         free arrays
    ****************************************************/

    for (i=0; i<n2; i++){
        free(S1[i]);
        free(S2[i]);
    }
    free(S1);
    free(S2);

}



#pragma optimization_level 1
void Hamiltonian_Band_Wannier(double ****RH, dcomplex **H, int Natom, int *FNAN, int **natn, int **ncn, int **atv_ijk, int *Total_NumOrbs, int *MP, double k1, double k2, double k3)
{
    int i,j,wanA,wanB,tnoA,tnoB,Anum,Bnum,NUM,GA_AN,LB_AN,GB_AN;
    int l1,l2,l3,Rn,n2;
    double **H1,**H2;
    double kRn,si,co,h;

    Anum = 1;
    for (i=0; i<Natom; i++){
        MP[i] = Anum;
        Anum += Total_NumOrbs[i];
    }
    NUM = Anum - 1;

    /****************************************************
                         Allocation
    ****************************************************/

    n2 = NUM + 2;

    H1 = (double**)malloc(sizeof(double*)*n2);
    for (i=0; i<n2; i++){
        H1[i] = (double*)malloc(sizeof(double)*n2);
    }

    H2 = (double**)malloc(sizeof(double*)*n2);
    for (i=0; i<n2; i++){
        H2[i] = (double*)malloc(sizeof(double)*n2);
    }

    /****************************************************
                        set Hamiltonian
    ****************************************************/

    H[0][0].r = 2.0*NUM;
    for (i=1; i<=NUM; i++){
        for (j=1; j<=NUM; j++){
            H1[i][j] = 0.0;
            H2[i][j] = 0.0;
        }
    }

    for (GA_AN=0; GA_AN<Natom; GA_AN++){
        tnoA = Total_NumOrbs[GA_AN];
        Anum = MP[GA_AN];

        for (LB_AN=0; LB_AN<=FNAN[GA_AN]; LB_AN++){
            GB_AN = natn[GA_AN][LB_AN]-1;
            Rn = ncn[GA_AN][LB_AN];
            tnoB = Total_NumOrbs[GB_AN];

            l1 = atv_ijk[Rn][0];
            l2 = atv_ijk[Rn][1];
            l3 = atv_ijk[Rn][2];
            kRn = k1*(double)l1 + k2*(double)l2 + k3*(double)l3;

            si = sin(2.0*PI*kRn);
            co = cos(2.0*PI*kRn);
            Bnum = MP[GB_AN];
            for (i=0; i<tnoA; i++){
                for (j=0; j<tnoB; j++){
                    h = RH[GA_AN][LB_AN][i][j];
                    H1[Anum+i][Bnum+j] += h*co;
                    H2[Anum+i][Bnum+j] += h*si;
                }
            }
        }
    }

    for (i=1; i<=NUM; i++){
        for (j=1; j<=NUM; j++){
            H[i][j].r = H1[i][j];
            H[i][j].i = H2[i][j];
        }
    }

    for (i=0; i<n2; i++){
        free(H1[i]);
        free(H2[i]);
    }
    free(H1);
    free(H2);
}



#pragma optimization_level 1
void Hamiltonian_Band_NC_Wannier(double *****RH, double *****IH, dcomplex **H, int Natom, int *FNAN, int **natn, int **ncn, int **atv_ijk, int *Total_NumOrbs, int *MP, double k1, double k2, double k3)
{
    int i,j,k,wanA,wanB,tnoA,tnoB,Anum,Bnum;
    int NUM,GA_AN,LB_AN,GB_AN;
    int l1,l2,l3,Rn,n2;
    double **H11r,**H11i;
    double **H22r,**H22i;
    double **H12r,**H12i;
    double kRn,si,co,h;

    /* set MP */

    Anum = 1;
    for (i=0; i<Natom; i++){
        MP[i] = Anum;
        Anum += Total_NumOrbs[i];
    }
    NUM = Anum - 1;
    n2 = NUM + 2;

    /*******************************************
     allocation of H11r, H11i,
                    H22r, H22i,
                    H12r, H12i
    *******************************************/

    H11r = (double**)malloc(sizeof(double*)*n2);
    for (i=0; i<n2; i++){
        H11r[i] = (double*)malloc(sizeof(double)*n2);
    }

    H11i = (double**)malloc(sizeof(double*)*n2);
    for (i=0; i<n2; i++){
        H11i[i] = (double*)malloc(sizeof(double)*n2);
    }

    H22r = (double**)malloc(sizeof(double*)*n2);
    for (i=0; i<n2; i++){
        H22r[i] = (double*)malloc(sizeof(double)*n2);
    }

    H22i = (double**)malloc(sizeof(double*)*n2);
    for (i=0; i<n2; i++){
        H22i[i] = (double*)malloc(sizeof(double)*n2);
    }

    H12r = (double**)malloc(sizeof(double*)*n2);
    for (i=0; i<n2; i++){
        H12r[i] = (double*)malloc(sizeof(double)*n2);
    }

    H12i = (double**)malloc(sizeof(double*)*n2);
    for (i=0; i<n2; i++){
        H12i[i] = (double*)malloc(sizeof(double)*n2);
    }

    /****************************************************
                        set Hamiltonian
    ****************************************************/

    H[0][0].r = 2.0*NUM;
    for (i=1; i<=NUM; i++){
        for (j=1; j<=NUM; j++){
            H11r[i][j] = 0.0;
            H11i[i][j] = 0.0;
            H22r[i][j] = 0.0;
            H22i[i][j] = 0.0;
            H12r[i][j] = 0.0;
            H12i[i][j] = 0.0;
        }
    }

    for (GA_AN=0; GA_AN<Natom; GA_AN++){
        tnoA = Total_NumOrbs[GA_AN];
        Anum = MP[GA_AN];

        for (LB_AN=0; LB_AN<=FNAN[GA_AN]; LB_AN++){
            GB_AN = natn[GA_AN][LB_AN]-1;
            Rn = ncn[GA_AN][LB_AN];
            tnoB = Total_NumOrbs[GB_AN];

            l1 = atv_ijk[Rn][0];
            l2 = atv_ijk[Rn][1];
            l3 = atv_ijk[Rn][2];
            kRn = k1*(double)l1 + k2*(double)l2 + k3*(double)l3;

            si = sin(2.0*PI*kRn);
            co = cos(2.0*PI*kRn);
            Bnum = MP[GB_AN];

            for (i=0; i<tnoA; i++){
                for (j=0; j<tnoB; j++){
                    H11r[Anum+i][Bnum+j] += co*RH[0][GA_AN][LB_AN][i][j] -  si*IH[0][GA_AN][LB_AN][i][j];
                    H11i[Anum+i][Bnum+j] += si*RH[0][GA_AN][LB_AN][i][j] +  co*IH[0][GA_AN][LB_AN][i][j];
                    H22r[Anum+i][Bnum+j] += co*RH[1][GA_AN][LB_AN][i][j] -  si*IH[1][GA_AN][LB_AN][i][j];
                    H22i[Anum+i][Bnum+j] += si*RH[1][GA_AN][LB_AN][i][j] +  co*IH[1][GA_AN][LB_AN][i][j];
                    H12r[Anum+i][Bnum+j] += co*RH[2][GA_AN][LB_AN][i][j] - si*(RH[3][GA_AN][LB_AN][i][j] + IH[2][GA_AN][LB_AN][i][j]);
                    H12i[Anum+i][Bnum+j] += si*RH[2][GA_AN][LB_AN][i][j] + co*(RH[3][GA_AN][LB_AN][i][j] + IH[2][GA_AN][LB_AN][i][j]);
                }
            }
        }
    }

    /******************************************************
        the full complex matrix of H
    ******************************************************/

    for (i=1; i<=NUM; i++){
        for (j=1; j<=NUM; j++){
            H[i    ][j    ].r =  H11r[i][j];
            H[i    ][j    ].i =  H11i[i][j];
            H[i+NUM][j+NUM].r =  H22r[i][j];
            H[i+NUM][j+NUM].i =  H22i[i][j];
            H[i    ][j+NUM].r =  H12r[i][j];
            H[i    ][j+NUM].i =  H12i[i][j];
            H[j+NUM][i    ].r =  H[i][j+NUM].r;
            H[j+NUM][i    ].i = -H[i][j+NUM].i;
        }
    }

    /****************************************************
                         free arrays
    ****************************************************/

    for (i=0; i<n2; i++){
        free(H11r[i]);
    }
    free(H11r);

    for (i=0; i<n2; i++){
        free(H11i[i]);
    }
    free(H11i);

    for (i=0; i<n2; i++){
        free(H22r[i]);
    }
    free(H22r);

    for (i=0; i<n2; i++){
        free(H22i[i]);
    }
    free(H22i);

    for (i=0; i<n2; i++){
        free(H12r[i]);
    }
    free(H12r);

    for (i=0; i<n2; i++){
        free(H12i[i]);
    }
    free(H12i);
}


#pragma optimization_level 1
void EigenState_k(
    double k1, double k2, double k3, int *MP, int spinsize, int SpinP_switch, 
    int Natom, int *FNAN, int **natn, int **ncn, int **atv_ijk, int *Total_NumOrbs,
    int fsize, int fsize2, int fsize3,  dcomplex ***Wk1, double **EigenVal1, 
    double ****OLP, double *****Hks, double *****iHks)
{
    int spin;
    int ik,i1,j1,i,j,l,k;
    int ct_AN,h_AN,mu1,mu2;
    int Anum,Bnum,tnoA,tnoB;
    int mn,jj1,ii1,m;
    int Rnh,Gh_AN,l1,l2,l3;
    double kx,ky,kz;
    
    double sumr,sumi;
    double tmp1r,tmp1i;
    double tmp2r,tmp2i;
    double tmp3r,tmp3i;
    double si,co,kRn,k1da,k1db,k1dc;

    double *ko,*M1;
    dcomplex **S,**H,**C; 
    
    double OLP_eigen_cut = 1.0e-12;
    dcomplex Ctmp1,Ctmp2;  
    int numprocs,myid,ID,ID1;


    /****************************************************
                     allocation of arrays:
    ****************************************************/

    ko = (double*)malloc(sizeof(double)*fsize3);
    M1 = (double*)malloc(sizeof(double)*fsize3);

    S = (dcomplex**)malloc(sizeof(dcomplex*)*fsize3);
    for (i=0; i<fsize3; i++){
        S[i] = (dcomplex*)malloc(sizeof(dcomplex)*fsize3);
        for (j=0; j<fsize3; j++){
            S[i][j].r = 0.0; 
            S[i][j].i = 0.0;
        }
    }

    H = (dcomplex**)malloc(sizeof(dcomplex*)*fsize3);
    for (i=0; i<fsize3; i++){
        H[i] = (dcomplex*)malloc(sizeof(dcomplex)*fsize3);
        for (j=0; j<fsize3; j++){ 
            H[i][j].r = 0.0;
            H[i][j].i = 0.0;
        }
    }

    C = (dcomplex**)malloc(sizeof(dcomplex*)*fsize3);
    for (i=0; i<fsize3; i++){
        C[i] = (dcomplex*)malloc(sizeof(dcomplex)*fsize3);
        for (j=0; j<fsize3; j++){
            C[i][j].r = 0.0;
            C[i][j].i = 0.0;
        }
    }
    

    Overlap_Band_Wannier(OLP,S,Natom,FNAN,natn,ncn,atv_ijk,Total_NumOrbs,MP,k1,k2,k3);
    EigenBand_lapack(S,ko,fsize,fsize,1);

    for (l=1; l<=fsize; l++){
        if (ko[l]<OLP_eigen_cut){
            printf("found an overcomplete basis set\n");
            exit(0); 
        }
    } 

    for (l=1; l<=fsize; l++){
        M1[l] = 1.0/sqrt(ko[l]);
    }

    for (i1=1; i1<=fsize; i1++){
        for (j1=1; j1<=fsize; j1++){
            S[i1][j1].r = S[i1][j1].r*M1[j1];
            S[i1][j1].i = S[i1][j1].i*M1[j1];
        } 
    } 

    for (spin=0; spin<spinsize; spin++){

        for (i1=1; i1<=fsize; i1++){
            for (j1=i1+1; j1<=fsize; j1++){
                Ctmp1 = S[i1][j1];
                Ctmp2 = S[j1][i1];
                S[i1][j1] = Ctmp2;
                S[j1][i1] = Ctmp1;
            }
        }

        if (SpinP_switch==0 || SpinP_switch==1){

            Hamiltonian_Band_Wannier(Hks[spin],H,Natom,FNAN,natn,ncn,atv_ijk,Total_NumOrbs,MP,k1,k2,k3);

            for (j1=1; j1<=fsize; j1++){
                for (i1=1; i1<=fsize; i1++){
                    sumr = 0.0;
                    sumi = 0.0;
                    for (l=1; l<=fsize; l++){
                        sumr += H[i1][l].r*S[j1][l].r - H[i1][l].i*S[j1][l].i;
                        sumi += H[i1][l].r*S[j1][l].i + H[i1][l].i*S[j1][l].r;
                    }
                    C[j1][i1].r = sumr;
                    C[j1][i1].i = sumi;
                }
            }     


            for (i1=1; i1<=fsize; i1++){
                for (j1=1; j1<=fsize; j1++){
                    sumr = 0.0;
                    sumi = 0.0;
                    for (l=1; l<=fsize; l++){
                        sumr +=  S[i1][l].r*C[j1][l].r + S[i1][l].i*C[j1][l].i;
                        sumi +=  S[i1][l].r*C[j1][l].i - S[i1][l].i*C[j1][l].r;
                    }
                    H[i1][j1].r = sumr;
                    H[i1][j1].i = sumi;
                }
            } 


            for (i1=1; i1<=fsize; i1++){
                for (j1=1; j1<=fsize; j1++){
                    C[i1][j1] = H[i1][j1];
                }
            }


            EigenBand_lapack(C,EigenVal1[spin],fsize,fsize,1);



            for (i1=1; i1<=fsize; i1++){
                for (j1=i1+1; j1<=fsize; j1++){
                    Ctmp1 = S[i1][j1];
                    Ctmp2 = S[j1][i1];
                    S[i1][j1] = Ctmp2;
                    S[j1][i1] = Ctmp1;
                }
            }

            for (i1=1; i1<=fsize; i1++){
                for (j1=i1+1; j1<=fsize; j1++){
                    Ctmp1 = C[i1][j1];
                    Ctmp2 = C[j1][i1];
                    C[i1][j1] = Ctmp2;
                    C[j1][i1] = Ctmp1;
                }
            }

            /* calculate wave functions */

            for (i1=1; i1<=fsize; i1++){
                for (j1=1; j1<=fsize; j1++){
                    sumr = 0.0;
                    sumi = 0.0;
                    for (l=1; l<=fsize; l++){
                        sumr +=  S[i1][l].r*C[j1][l].r - S[i1][l].i*C[j1][l].i;
                        sumi +=  S[i1][l].r*C[j1][l].i + S[i1][l].i*C[j1][l].r;
                    }
                    Wk1[spin][j1][i1].r = sumr;
                    Wk1[spin][j1][i1].i = sumi;
                }
            }
        }else if (SpinP_switch==3){

            Hamiltonian_Band_NC_Wannier(Hks,iHks,H,Natom,FNAN,natn,ncn,atv_ijk,Total_NumOrbs,MP,k1,k2,k3);

            for (j1=1; j1<=fsize; j1++){
                for (i1=1; i1<=fsize2; i1++){
                    for (m=0; m<=1; m++){
                        sumr = 0.0;
                        sumi = 0.0;
                        mn = m*fsize;
                        for (l=1; l<=fsize; l++){
                            sumr += H[i1][l+mn].r*S[j1][l].r - H[i1][l+mn].i*S[j1][l].i;
                            sumi += H[i1][l+mn].r*S[j1][l].i + H[i1][l+mn].i*S[j1][l].r;
                        }
                        jj1 = 2*j1 - 1 + m;
                        C[jj1][i1].r = sumr;
                        C[jj1][i1].i = sumi;
                    }
                }
            }     

            for (i1=1; i1<=fsize; i1++){
                for (m=0; m<=1; m++){
                    ii1 = 2*i1 - 1 + m;
                    for (j1=1; j1<=fsize2; j1++){
                        sumr = 0.0;
                        sumi = 0.0;
                        mn = m*fsize;
                        for (l=1; l<=fsize; l++){
                            sumr +=  S[i1][l].r*C[j1][l+mn].r + S[i1][l].i*C[j1][l+mn].i;
                            sumi +=  S[i1][l].r*C[j1][l+mn].i - S[i1][l].i*C[j1][l+mn].r;
                        }
                        H[ii1][j1].r = sumr;
                        H[ii1][j1].i = sumi;
                    }
                }
            }


            EigenBand_lapack(H,EigenVal1[0],fsize2,fsize2,1);


            for (i1=1; i1<=fsize; i1++){
                for (j1=i1+1; j1<=fsize; j1++){
                    Ctmp1 = S[i1][j1];
                    Ctmp2 = S[j1][i1];
                    S[i1][j1] = Ctmp2;
                    S[j1][i1] = Ctmp1;
                }
            }

            for (i1=1; i1<=fsize2; i1++){
                for (j1=1; j1<=fsize2; j1++){
                    C[i1][j1].r = 0.0;
                    C[i1][j1].i = 0.0;
                }
            }

            for (m=0; m<=1; m++){
                for (i1=1; i1<=fsize; i1++){
                    for (j1=1; j1<=fsize2; j1++){
                        sumr = 0.0; 
                        sumi = 0.0;
                        for (l=1; l<=fsize; l++){
                            sumr +=  S[i1][l].r*H[2*l-1+m][j1].r - S[i1][l].i*H[2*l-1+m][j1].i;
                            sumi +=  S[i1][l].r*H[2*l-1+m][j1].i + S[i1][l].i*H[2*l-1+m][j1].r;
                        } 
                        Wk1[0][j1][i1+m*fsize].r = sumr;
                        Wk1[0][j1][i1+m*fsize].i = sumi;
                    }
                }
            }
        }
    }

    
    for (i=0; i<fsize3; i++){
        free(C[i]);
    }
    free(C);
    
    for (i=0; i<fsize3; i++){
        free(H[i]); 
    }
    free(H);
    
    for (i=0; i<fsize3; i++){
        free(S[i]); 
    }    
    free(S);

    free(M1); 
    free(ko);
}