#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
#include <math.h>
#include "Generate_MLWF.h"


void Write_win(
    char *filename, char **Atoms_Symbol, int Natom, double **Gxyz, double **tv, double **rtv, int BANDNUM, int Wannier_Func_Num, 
    double Wannier_Outer_Window_Bottom, double Wannier_Outer_Window_Top, double Wannier_Inner_Window_Bottom, double Wannier_Inner_Window_Top,
    int knum_i, int knum_j, int knum_k, double **kg)
{
    int i, k, Gc_AN, itmp;
    double Cxyz[3], Cell_Gxyz[3];
    int kpt_num;
    char fname[300];    
    FILE *fp;

    kpt_num = knum_i*knum_j*knum_k;

    sprintf(fname,"%s.win",filename);
    if((fp=fopen(fname,"wt"))==NULL){
        printf("Error in opening %s for input file for wannier90.\n",fname);
        exit(0);
    }

    fprintf(fp, "num_bands %d\n", BANDNUM);
    fprintf(fp, "num_wann %d\n", Wannier_Func_Num);
    fprintf(fp, "\n");
    fprintf(fp, "dis_win_min =%18.14f\n",Wannier_Outer_Window_Bottom); 
    fprintf(fp, "dis_win_max =%18.14f\n",Wannier_Outer_Window_Top); 
    fprintf(fp, "dis_froz_min =%18.14f\n",Wannier_Inner_Window_Bottom); 
    fprintf(fp, "dis_froz_max =%18.14f\n",Wannier_Inner_Window_Top); 
    fprintf(fp, "\n");
    fprintf(fp, "dis_num_iter 0\n");
    fprintf(fp, "num_iter 0\n");
    fprintf(fp, "\n");
    fprintf(fp, "begin unit_cell_cart\n");
    fprintf(fp, "Ang\n");
    for (i=0; i<=2; i++){
        fprintf(fp,"%lf %lf %lf\n",tv[i][0]*BohrR_Wannier,tv[i][1]*BohrR_Wannier,tv[i][2]*BohrR_Wannier);
    }
    fprintf(fp, "end unit_cell_cart\n");
    fprintf(fp, "\n");
    fprintf(fp, "begin atoms_frac\n");
    for (Gc_AN=0; Gc_AN<Natom; Gc_AN++){

        Cxyz[0] = Gxyz[Gc_AN][0];
        Cxyz[1] = Gxyz[Gc_AN][1];
        Cxyz[2] = Gxyz[Gc_AN][2];

        Cell_Gxyz[0] = (Cxyz[0]*rtv[0][0]+Cxyz[1]*rtv[0][1]+Cxyz[2]*rtv[0][2])*0.5/PI;
        Cell_Gxyz[1] = (Cxyz[0]*rtv[1][0]+Cxyz[1]*rtv[1][1]+Cxyz[2]*rtv[1][2])*0.5/PI;
        Cell_Gxyz[2] = (Cxyz[0]*rtv[2][0]+Cxyz[1]*rtv[2][1]+Cxyz[2]*rtv[2][2])*0.5/PI;

        for (i=0; i<=2; i++){

            itmp = (int)Cell_Gxyz[i]; 

            if (1.0<Cell_Gxyz[i]){
                Cell_Gxyz[i] = fabs(Cell_Gxyz[i] - (double)itmp);
            }else if (Cell_Gxyz[i]<-1.0e-13){
                Cell_Gxyz[i] = fabs(Cell_Gxyz[i] + (double)(abs(itmp)+1));
            }
        }

        fprintf(fp,"%4s %18.14f %18.14f %18.14f\n", Atoms_Symbol[Gc_AN], Cell_Gxyz[0], Cell_Gxyz[1], Cell_Gxyz[2]);
    }
    fprintf(fp, "end atoms_frac\n");
    fprintf(fp, "\n");
    fprintf(fp, "mp_grid = %d %d %d\n",knum_i, knum_j, knum_k);
    fprintf(fp, "begin kpoints\n");
    for(k=0;k<kpt_num;k++){
        fprintf(fp, "%12.8f%12.8f%12.8f\n",kg[k][0],kg[k][1],kg[k][2]);
    }
    fprintf(fp, "End Kpoints\n");
}


void Write_eig(char *filename, int BANDNUM, int Nkpt, int spinsize, int fsize3, int ***Nk, double ***EigenValall, double ChemP)
{
    int spin, k, i;
    double dtmp = 10000.0;
    char fname[300];
    FILE *fp;

    sprintf(fname,"%s.eig",filename);

    if((fp=fopen(fname,"wt"))==NULL){
        printf("Error in opening %s for writing eigenvalues.\n",fname);
        exit(0);
    }

    for(spin=0; spin<spinsize; spin++){
        for(k=0; k<Nkpt; k++){
            for(i=1; i<BANDNUM+1; i++){
                if ( (i+Nk[spin][k][0])<fsize3 && (i+Nk[spin][k][0])<=Nk[spin][k][1]){
                    fprintf(fp,"%5d%5d%18.12f\n",i,k+1,(EigenValall[k][spin][i+Nk[spin][k][0]]-ChemP)*eV2Hartree);
                }else {
                    fprintf(fp,"%5d%5d  %18.10f\n",i,k+1,dtmp);
                }
            }
        }
    }

    fclose(fp);
}


void Write_Amnk(char *filename, int BANDNUM, int Nkpt, int WANNUM, int spinsize, dcomplex ****Amnk)
{
    int spin, k, nindx, mu1;
    char fname[300];
    FILE *fp;

    sprintf(fname,"%s.amn",filename);

    if((fp=fopen(fname,"wt"))==NULL){
        printf("******************************************************************\n");
        printf("* Error in opening file for Amn(k).\n");
        printf("******************************************************************\n");
    }else{
        printf(" ... ... Writting Amn(k) matrix into file.\n\n");
        fprintf(fp,"Amn. Fist line BANDNUM, KPTNUM, WANNUM, spinsize. Next is m n k and elements.Spin is the most outer loop.\n");
        fprintf(fp,"%13d%13d%13d%13d\n",BANDNUM,Nkpt,WANNUM,spinsize);
    }


    for(spin=0;spin<spinsize; spin++){
        for(k=0;k<Nkpt;k++){
            for(nindx=0;nindx<WANNUM;nindx++){
                for(mu1=0;mu1<BANDNUM;mu1++){
                    fprintf(fp,"%5d%5d%5d%18.12f%18.12f\n",mu1+1,nindx+1,k+1,Amnk[spin][k][mu1][nindx].r,Amnk[spin][k][mu1][nindx].i);
                }
            }
        }
    }

    fclose(fp);
}


void Write_Mmnkb(
    char *filename, 
    int BANDNUM, int Nkpt, int tot_bvector, int spinsize,
    double **kg, double **frac_bv, int **kplusb, dcomplex *****Mmnkb_zero)
{
    int spin, k, bindx, kk, i, j;
    double b[3], ktp[3];
    char fname[300];
    FILE *fp;

    sprintf(fname,"%s.mmn",filename);

    if((fp=fopen(fname,"wt"))==NULL){
        printf("******************************************************************\n");
        printf("* Error in opening file %s for writing Mmn(k,b).\n",fname);
        printf("******************************************************************\n");
        exit(0);     
    }else{
        printf(" ... ... Writting Mmn_zero(k,b) matrix into file.\n\n");
        fprintf(fp,"Mmn_zero(k,b). band_num, kpt_num, bvector num, spinsize\n");
        fprintf(fp,"%13d%13d%13d%13d\n",BANDNUM,Nkpt,tot_bvector,spinsize);    
    }

    for(spin=0;spin<spinsize; spin++){
        for(k=0;k<Nkpt;k++){
            for(bindx=0;bindx<tot_bvector;bindx++){ 

                kk=kplusb[k][bindx];
                ktp[0]=kg[k][0]+frac_bv[bindx][0];  /* Equivalent k point of k+b */
                ktp[1]=kg[k][1]+frac_bv[bindx][1];
                ktp[2]=kg[k][2]+frac_bv[bindx][2];

                for(i=0;i<=2;i++){
                    b[i]=ktp[i];
                    if(ktp[i]>=1.0){
                        b[i]=ktp[i]-1.0;
                    }
                    if(ktp[i]<0.0){
                        b[i]=ktp[i]+1.0;
                    }
                }

                fprintf(fp,"%5d%5d%5d%5d%5d\n", k+1,kk+1, (int)(ktp[0]-b[0]), (int)(ktp[1]-b[1]), (int)(ktp[2]-b[2])); 

                for(i=1;i<BANDNUM+1;i++){
                    for(j=1;j<BANDNUM+1;j++){
                        fprintf(fp,"%18.12f%18.12f\n", Mmnkb_zero[k][bindx][spin][j][i].r, Mmnkb_zero[k][bindx][spin][j][i].i);
                    }
                }
            }
        }
    }
    fclose(fp);
}
