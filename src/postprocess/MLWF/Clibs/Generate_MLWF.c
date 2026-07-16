#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include "Generate_MLWF.h"


void LCPAO2Wannier(
    char *filename,
    char **Atoms_Symbol, int SpinP_switch, int Natom, double **Gxyz, int *FNAN, int **natn, int **ncn, 
    int **atv_ijk, int *Total_NumOrbs, double **Latvecs, double **Recvecs,
    int knum_i, int knum_j, int knum_k, int Nkpt, int WANNUM,
    int tot_bvector, double **bvector, double *wb, int **kplusb, double **kg, double **frac_bv,
    int MLWF_Num_Kinds_Projectors, int **MLWF_NumL_Pro,
    int *MLWF_Num_Pro, int **MLWF_Select_Matrix, double ***MLWF_Projector_Hybridize_Matrix, double ****MLWF_RotMat_for_Real_Func,
    double MLWF_Outer_Window_Bottom, double MLWF_Outer_Window_Top, double MLWF_Inner_Window_Bottom, double MLWF_Inner_Window_Top,
    int MLWF_Dis_SCF_Max_Steps, double MLWF_Dis_Conv_Criterion, double MLWF_Dis_Mixing_Para,
    int MLWF_Min_Secant_Steps, double MLWF_Min_StepLength, int MLWF_Min_Scheme, int MLWF_Minimizing_Max_Steps, double MLWF_Min_Conv_Criterion, double MLWF_Min_Secant_StepLength,
    double **Wannier_Guide,
    int *FNAN_WP, int **natn_WP, int **ncn_WP,
    double ****OLP_WP,
    double ****OLP, double *****Hks, double *****iHks, dcomplex *****OLPe, 
    double ChemP)
{
    
    int spinsize, fsize, fsize2, fsize3;
    int spin, i, j, itmp1, k, ki, kj, itmp2, nindx, bindx, itmp3, m, n, kk, m1, n1, l;
    int odloop, odloop_num;
    double k1[3], k2[3], bk[3], dk[3], Sop[2];
    double norm, wbtot, sum, sumr;
    double oemin, oemax, iemin, iemax;
    double tmpr, tmpi;
    int OMPID, Nthrds, Nprocs;
    int BANDNUM, MkNUM;
    int *MP;

    dcomplex ****Wkall;
    dcomplex *****Mmnkb_zero;
    dcomplex *****Mmnkb_dis;
    dcomplex ****M_zero,****Mmnkb; 
    dcomplex ****Amnk, ****Uk, ****Utilde;
    dcomplex ***deltaW,***Ukmatrix, **tmpM;
    dcomplex ***Gmnk, *rmnk, *dmnk, *rimnk;
    dcomplex ***csheet;
    double ***sheet;
    double **rguide;
    double ***EigenValall;
    int ***Nk; 
    int ***innerNk;
    double ***eigen;

    double **wann_center;
    double *wann_r2;
    double sumr2;
    double omega, omega_I, omega_OD, omega_D;
    double omega_prev, omega_I_prev, omega_OD_prev, omega_D_prev;
    double omega_next, omega_I_next, omega_OD_next, omega_D_next;
    double delta_new, delta_0, delta_old, delta_mid, delta_d;
    double yita,  yita_prev,  beta;
    double epcg, epsec, sigma_0; 

    double alpha,conv_total_omega;
    int step;
    int sec_step, sec_max;

    int searching_scheme;
    int max_steps = MLWF_Minimizing_Max_Steps;



    MP = (int*)malloc(sizeof(int)*Natom);
  
    fsize = 1;
    for (i=0; i<Natom; i++){
        MP[i] = fsize;
        fsize += Total_NumOrbs[i];
    }
    fsize--;
    
    if      (SpinP_switch==0){ spinsize=1; fsize2=fsize;  fsize3=fsize+2;}
    else if (SpinP_switch==1){ spinsize=2; fsize2=fsize;  fsize3=fsize+2;}
    else if (SpinP_switch==3){ spinsize=1; fsize2=2*fsize;fsize3=2*fsize+2;}

    wbtot = 0.0;
    for (bindx=0; bindx<tot_bvector; bindx++){
        wbtot=wbtot+wb[bindx];
    }


    /* for entangled bands */
    Nk = (int***)malloc(sizeof(int**)*spinsize);
    for(spin=0;spin<spinsize;spin++){
        Nk[spin]=(int**)malloc(sizeof(int*)*Nkpt);
        for(k=0;k<Nkpt;k++){
            Nk[spin][k]=(int*)malloc(sizeof(int)*2);
        } 
    }

    innerNk = (int***)malloc(sizeof(int**)*spinsize);
    for(spin=0;spin<spinsize;spin++){
        innerNk[spin]=(int**)malloc(sizeof(int*)*Nkpt);
        for(k=0;k<Nkpt;k++){
            innerNk[spin][k]=(int*)malloc(sizeof(int)*2);
        }
    }
    
    Wkall = (dcomplex****)malloc(sizeof(dcomplex***)*Nkpt); 
    for(k=0;k<Nkpt;k++){
        Wkall[k] = (dcomplex***)malloc(sizeof(dcomplex**)*spinsize); 
        for (spin=0; spin<spinsize; spin++){
            Wkall[k][spin] = (dcomplex**)malloc(sizeof(dcomplex*)*(fsize3)); 
            for (i=0; i<fsize3; i++){
                Wkall[k][spin][i] = (dcomplex*)malloc(sizeof(dcomplex)*fsize3); 
                for (j=0; j<fsize3; j++){ 
                    Wkall[k][spin][i][j].r=0.0;
                    Wkall[k][spin][i][j].i=0.0;
                }
            } 
        }
    }
    

    EigenValall = (double***)malloc(sizeof(double**)*Nkpt);
    for(k=0;k<Nkpt;k++){
        EigenValall[k] = (double**)malloc(sizeof(double*)*spinsize);
        for (spin=0; spin<spinsize; spin++){
            EigenValall[k][spin] = (double*)malloc(sizeof(double)*(fsize3));
            for (j=0; j<fsize3; j++){
                EigenValall[k][spin][j] = 1.0e9;
            }
        }
    }

    /*
    #pragma omp parallel shared(Nkpt,kg,MP,spinsize,SpinP_switch,fsize,fsize2,fsize3,Wkall,EigenValall,OLP,Hks,iHks) private(k,OMPID,Nthrds)
    {
        OMPID = omp_get_thread_num();
        Nthrds = omp_get_num_threads();
        Nprocs = omp_get_num_procs();

        for( k=OMPID; k<Nkpt; k+=Nthrds ){
            EigenState_k(
                kg[k][0], kg[k][1], kg[k][2], 
                MP, spinsize, SpinP_switch, 
                Natom, FNAN, natn, ncn, atv_ijk, Total_NumOrbs, 
                fsize, fsize2, fsize3, 
                Wkall[k], EigenValall[k], OLP, Hks, iHks );
        }
    }*/


    for( k=0; k<Nkpt; k+=1 ){
        EigenState_k(
            kg[k][0], kg[k][1], kg[k][2], 
            MP, spinsize, SpinP_switch, 
            Natom, FNAN, natn, ncn, atv_ijk, Total_NumOrbs, 
            fsize, fsize2, fsize3, 
            Wkall[k], EigenValall[k], OLP, Hks, iHks );
    }



    oemin=MLWF_Outer_Window_Bottom/eV2Hartree+ChemP;
    oemax=MLWF_Outer_Window_Top/eV2Hartree+ChemP;
    iemin=MLWF_Inner_Window_Bottom/eV2Hartree+ChemP;
    iemax=MLWF_Inner_Window_Top/eV2Hartree+ChemP;

    BANDNUM = 0; 
    MkNUM = 0;
    
    printf("Selected bands within the Outer and Inner Windows at each k point:\n");
    printf("spin| kpt | Outer Window (Nk) |  Inner Window (Mk)  |\n"); 


    for(k=0; k<Nkpt; k++){

        for(spin=0; spin<spinsize; spin++){
            
            printf(" %1d  |%5d|",spin,k+1);

            Nk[spin][k][0]      = 0;
            Nk[spin][k][1]      = 0;
            innerNk[spin][k][0] = 0;
            innerNk[spin][k][1] = 0;

            for(i=1; i<fsize3-1; i++){

                if (EigenValall[k][spin][i]<oemin){
                    Nk[spin][k][0] = i;
                }
                if (EigenValall[k][spin][i]<oemax){
                    Nk[spin][k][1] = i;
                }
                if (EigenValall[k][spin][i]<iemin){
                    innerNk[spin][k][0] = i;
                }
                if (EigenValall[k][spin][i]<iemax){
                    innerNk[spin][k][1] = i;
                } 
            } 

            if( BANDNUM<(Nk[spin][k][1]-Nk[spin][k][0]) ){
                BANDNUM = Nk[spin][k][1] - Nk[spin][k][0];
            }

            if( MkNUM<(innerNk[spin][k][1]-innerNk[spin][k][0]) ){
                MkNUM = innerNk[spin][k][1] - innerNk[spin][k][0];
            }

            printf("  (%3d ,%3d]  %3d  |  (%3d ,%3d]   %3d   |\n",
                    Nk[spin][k][0],Nk[spin][k][1],
                    Nk[spin][k][1]-Nk[spin][k][0],
                    innerNk[spin][k][0],
                    innerNk[spin][k][1],
                    innerNk[spin][k][1]-innerNk[spin][k][0]);
            fflush(0);

            if (WANNUM>(Nk[spin][k][1]-Nk[spin][k][0])){
                printf("**************************ERROR**************************\n");
                printf("* Bands number within OUTER window [%10.5f, %10.5f]     *\n",oemin,oemax);
                printf("* is less than Wannier Function number %i.              *\n",WANNUM);
                printf("* Please check them and try again.                      *\n");
                printf("**************************ERROR**************************\n");      
                exit(0); 
            }

            if (WANNUM<(innerNk[spin][k][1]-innerNk[spin][k][0])){
                printf("**************************ERROR**************************\n");
                printf("* Bands number within INNER window [%10.5f, %10.5f]     *\n",iemin,iemax);
                printf("* is larger than wannier function number %i.            *\n",WANNUM);
                printf("* Please check them and try again.\n");
                printf("**************************ERROR**************************\n");
                exit(0);
            }
        }
    }

    Write_eig(filename, BANDNUM, Nkpt, spinsize, fsize3, Nk, EigenValall, ChemP);
    Write_win(filename, Atoms_Symbol, Natom, Gxyz, Latvecs, Recvecs, BANDNUM, WANNUM, 
              MLWF_Outer_Window_Bottom, MLWF_Outer_Window_Top, MLWF_Inner_Window_Bottom, MLWF_Inner_Window_Top,
              knum_i, knum_j, knum_k, kg);



    eigen = (double***)malloc(sizeof(double**)*spinsize);

    for (spin=0;spin<spinsize;spin++){
        eigen[spin] = (double**)malloc(sizeof(double*)*Nkpt);
        for (k=0; k<Nkpt; k++){
            eigen[spin][k] = (double*)malloc(sizeof(double)*BANDNUM);
            for (j=0; j<BANDNUM; j++){ 
                eigen[spin][k][j] = 0.0;
            }
        }
    }

    for (spin=0;spin<spinsize;spin++){
        for (k=0;k<Nkpt;k++){
            for (i=1; i<BANDNUM+1; i++){
                eigen[spin][k][i-1] = EigenValall[k][spin][i+Nk[spin][k][0]];
            }
        }
    }


    for(k=0; k<Nkpt; k++){
        for (spin=0; spin<spinsize; spin++){
            free(EigenValall[k][spin]);
        }
        free(EigenValall[k]);
    }
    free(EigenValall);



    Mmnkb_zero = (dcomplex*****)malloc(sizeof(dcomplex****)*Nkpt);
    for(k=0;k<Nkpt;k++){
        Mmnkb_zero[k] = (dcomplex****)malloc(sizeof(dcomplex***)*tot_bvector);
        for(bindx=0;bindx<tot_bvector;bindx++){
            Mmnkb_zero[k][bindx] = (dcomplex***)malloc(sizeof(dcomplex**)*spinsize);
            for (spin=0; spin<spinsize; spin++){
                Mmnkb_zero[k][bindx][spin] = (dcomplex**)malloc(sizeof(dcomplex*)*(BANDNUM+1));
                for (i=0; i<BANDNUM+1; i++){
                    Mmnkb_zero[k][bindx][spin][i] = (dcomplex*)malloc(sizeof(dcomplex)*(BANDNUM+1));
                    for (j=0; j<BANDNUM+1; j++){ 
                        Mmnkb_zero[k][bindx][spin][i][j].r=0.0; 
                        Mmnkb_zero[k][bindx][spin][i][j].i=0.0;
                    }
                }
            }
        }
    }


    printf("\nComing to the overlap matrix calculating......\n");
    fflush(0);

    odloop_num = spinsize*Nkpt*tot_bvector*BANDNUM*BANDNUM;

    #pragma omp parallel shared(odloop_num,Nkpt,tot_bvector,BANDNUM,kplusb,Nk,fsize3,kg,frac_bv,SpinP_switch,MP,Recvecs,OLPe,Wkall,fsize,Mmnkb_zero) private(odloop,spin,itmp1,itmp2,itmp3,k,bindx,m,n,kk,m1,n1,k1,k2,bk,dk,Sop,norm,OMPID,Nthrds,Nprocs)
    {

        OMPID = omp_get_thread_num();
        Nthrds = omp_get_num_threads();
        Nprocs = omp_get_num_procs();

        for( odloop=OMPID; odloop<odloop_num; odloop+=Nthrds ){

            spin = odloop/(Nkpt*tot_bvector*BANDNUM*BANDNUM);
            itmp1 = odloop - spin*(Nkpt*tot_bvector*BANDNUM*BANDNUM);
            k = itmp1/(tot_bvector*BANDNUM*BANDNUM);
            itmp2 = itmp1 - k*(tot_bvector*BANDNUM*BANDNUM);
            bindx = itmp2/(BANDNUM*BANDNUM);
            itmp3 = itmp2 - bindx*(BANDNUM*BANDNUM);
            m = itmp3/BANDNUM;
            n = itmp3 - m*BANDNUM;

            kk = kplusb[k][bindx];

            m1 = m + 1 + Nk[spin][k][0];
            n1 = n + 1 + Nk[spin][kk][0];      

            if (m1<(fsize3-1) && n1<(fsize3-1)){

                k1[0] = kg[k][0];
                k1[1] = kg[k][1];
                k1[2] = kg[k][2];

                k2[0] = k1[0] + frac_bv[bindx][0];
                k2[1] = k1[1] + frac_bv[bindx][1];
                k2[2] = k1[2] + frac_bv[bindx][2];

                bk[0] = frac_bv[bindx][0]; 
                bk[1] = frac_bv[bindx][1];
                bk[2] = frac_bv[bindx][2];

                dk[0] = bk[0]*Recvecs[0][0] + bk[1]*Recvecs[1][0] + bk[2]*Recvecs[2][0];
                dk[1] = bk[0]*Recvecs[0][1] + bk[1]*Recvecs[1][1] + bk[2]*Recvecs[2][1];
                dk[2] = bk[0]*Recvecs[0][2] + bk[1]*Recvecs[1][2] + bk[2]*Recvecs[2][2];

                Calc_Mmnkb_zero(
                    k1, k2, dk, bindx,
                    SpinP_switch, Natom, Gxyz, FNAN, natn, ncn, 
                    atv_ijk, Total_NumOrbs, MP, fsize,
                    Sop, OLPe, Wkall[k][spin][m1], Wkall[kk][spin][n1]);

                Mmnkb_zero[k][bindx][spin][m+1][n+1].r = Sop[0];
                Mmnkb_zero[k][bindx][spin][m+1][n+1].i = Sop[1];

                norm = sqrt( fabs(Sop[0]*Sop[0] + Sop[1]*Sop[1]) );

                if (norm>1.0){
                    printf("**********************WARNNING**********************\n");
                    printf("Attention! |Mmnkb=%10.5f|>1.0 at k=%i,b=%i,i=%i,j=%i\n",norm,k,bindx,m+1,n+1);
                    printf("**********************WARNNING**********************\n");
                }
            }else {
                Mmnkb_zero[k][bindx][spin][m+1][n+1].r = -99999.0;
                Mmnkb_zero[k][bindx][spin][m+1][n+1].i = -99999.0;
            }
        }
    }

    Write_Mmnkb(filename, BANDNUM, Nkpt, tot_bvector, spinsize, kg, frac_bv, kplusb, Mmnkb_zero);


    Amnk=(dcomplex****)malloc(sizeof(dcomplex***)*spinsize);
    for(spin=0;spin<spinsize;spin++){
        Amnk[spin]=(dcomplex***)malloc(sizeof(dcomplex**)*Nkpt);
        for(k=0;k<Nkpt;k++){
            Amnk[spin][k]=(dcomplex**)malloc(sizeof(dcomplex*)*BANDNUM);
            for(i=0;i<BANDNUM;i++){
                Amnk[spin][k][i]=(dcomplex*)malloc(sizeof(dcomplex)*WANNUM);
                for(j=0;j<WANNUM;j++){
                    Amnk[spin][k][i][j].r=0.0;
                    Amnk[spin][k][i][j].i=0.0;
                    if(j==i){
                        Amnk[spin][k][i][j].r=1.0;
                        Amnk[spin][k][i][j].i=0.0;
                    }
                }
            }
        }
    }

    Projection_Amatrix(
        MLWF_Num_Kinds_Projectors, MLWF_NumL_Pro,
        MLWF_Num_Pro, MLWF_Select_Matrix, MLWF_Projector_Hybridize_Matrix, MLWF_RotMat_for_Real_Func,
        FNAN_WP, natn_WP, ncn_WP, atv_ijk, Total_NumOrbs, 
        OLP_WP,
        Amnk, kg, spinsize, 
        fsize, SpinP_switch, 
        Nkpt, BANDNUM, WANNUM, 
        Wkall, MP, Nk);

    Write_Amnk(filename, BANDNUM, Nkpt, WANNUM, spinsize, Amnk);
 

    for(spin=0;spin<spinsize;spin++){
        for(k=0;k<Nkpt;k++){
            free(innerNk[spin][k]);
        }
        free(innerNk[spin]);
    }
    free(innerNk);

    for(spin=0;spin<spinsize;spin++){
        for(k=0;k<Nkpt;k++){
            free(Nk[spin][k]);
        }
        free(Nk[spin]);  
    }
    free(Nk);

    for(k=0;k<Nkpt;k++){
        for (spin=0; spin<spinsize; spin++){
            for (i=0; i<fsize3; i++){
                free(Wkall[k][spin][i]);
            }
            free(Wkall[k][spin]);
        }
        free(Wkall[k]);
    }
    free(Wkall);

    for(k=0;k<Nkpt;k++){
        for(bindx=0;bindx<tot_bvector;bindx++){
            for(spin=0; spin<spinsize; spin++){
                for (i=0; i<BANDNUM+1; i++){
                    free(Mmnkb_zero[k][bindx][spin][i] );
                }
                free(Mmnkb_zero[k][bindx][spin]);
            }
            free(Mmnkb_zero[k][bindx]);
        }
        free(Mmnkb_zero[k]);
    }
    free(Mmnkb_zero);

    for(spin=0;spin<spinsize;spin++){
        for (k=0; k<Nkpt; k++){
            free(eigen[spin][k]);
        }
        free(eigen[spin]); 
    }
    free(eigen);

    free(MP);
}