#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "Generate_MLWF.h"

#pragma optimization_level 1
void Calc_Mmnkb_zero(
    double k1[3], double k2[3], double dk[3], int bindx,
    int SpinP_switch, int Natom, double **Gxyz, int *FNAN, int **natn, int **ncn, 
    int **atv_ijk, int *Total_NumOrbs, int *MP, int fsize, double Sop[2], 
    dcomplex *****OLPe, dcomplex *Wk1, dcomplex *Wk2)
{
    int ct_AN,tnoA,Anum,h_AN;
    int i1,j1,l1,l2,l3,Rnh,Gh_AN,tnoB,Bnum;
    double si1,co1,si2,co2,kRn;
    double tmp1r,tmp1i,tmp2r1,tmp2i1,tmp2r2,tmp2i2;
    double tmp3r,tmp3i,sumr,sumi;


    if (SpinP_switch==0 || SpinP_switch==1){

        sumr = 0.0; 
        sumi = 0.0; 

        for (ct_AN=0; ct_AN<Natom; ct_AN++){ 

            tnoA = Total_NumOrbs[ct_AN];
            Anum = MP[ct_AN];

            for (h_AN=0; h_AN<=FNAN[ct_AN]; h_AN++){

                Rnh = ncn[ct_AN][h_AN];
                Gh_AN = natn[ct_AN][h_AN]-1;
                tnoB = Total_NumOrbs[Gh_AN];
                Bnum = MP[Gh_AN];

                l1 = atv_ijk[Rnh][0];
                l2 = atv_ijk[Rnh][1];
                l3 = atv_ijk[Rnh][2];

                kRn = 2.0*PI*(k2[0]*(double)l1 + k2[1]*(double)l2 + k2[2]*(double)l3);
                kRn -= (dk[0]*Gxyz[ct_AN][0] + dk[1]*Gxyz[ct_AN][1] + dk[2]*Gxyz[ct_AN][2]);

                si1 = sin(kRn);
                co1 = cos(kRn);

                kRn = -2.0*PI*(k1[0]*(double)l1 + k1[1]*(double)l2 + k1[2]*(double)l3);
                kRn -= (dk[0]*Gxyz[ct_AN][0] + dk[1]*Gxyz[ct_AN][1] + dk[2]*Gxyz[ct_AN][2]);

                si2 = sin(kRn);
                co2 = cos(kRn);

                for (i1=0; i1<tnoA; i1++){
                    for (j1=0; j1<tnoB; j1++){

                        tmp1r = Wk1[Anum+i1].r*Wk2[Bnum+j1].r + Wk1[Anum+i1].i*Wk2[Bnum+j1].i;
                        tmp1i = Wk1[Anum+i1].r*Wk2[Bnum+j1].i - Wk1[Anum+i1].i*Wk2[Bnum+j1].r;

                        tmp2r1 = co1*tmp1r - si1*tmp1i;  
                        tmp2i1 = co1*tmp1i + si1*tmp1r; 

                        tmp1r = Wk1[Bnum+j1].r*Wk2[Anum+i1].r + Wk1[Bnum+j1].i*Wk2[Anum+i1].i;
                        tmp1i = Wk1[Bnum+j1].r*Wk2[Anum+i1].i - Wk1[Bnum+j1].i*Wk2[Anum+i1].r;

                        tmp2r2 = co2*tmp1r - si2*tmp1i;
                        tmp2i2 = co2*tmp1i + si2*tmp1r;

                        tmp3r = OLPe[bindx][ct_AN][h_AN][i1][j1].r;
                        tmp3i = OLPe[bindx][ct_AN][h_AN][i1][j1].i;

                        sumr += 0.5*(tmp2r1+tmp2r2)*tmp3r - 0.5*(tmp2i1+tmp2i2)*tmp3i;
                        sumi += 0.5*(tmp2r1+tmp2r2)*tmp3i + 0.5*(tmp2i1+tmp2i2)*tmp3r;
                    }
                }
            }   
        }       

        Sop[0] = sumr;
        Sop[1] = sumi;

    }else if (SpinP_switch==3){

        sumr = 0.0; 
        sumi = 0.0; 

        for (ct_AN=0; ct_AN<Natom; ct_AN++){ 

            tnoA = Total_NumOrbs[ct_AN];
            Anum = MP[ct_AN];

            for (h_AN=0; h_AN<=FNAN[ct_AN]; h_AN++){

                Rnh = ncn[ct_AN][h_AN];
                Gh_AN = natn[ct_AN][h_AN]-1;
                tnoB = Total_NumOrbs[Gh_AN];
                Bnum = MP[Gh_AN];

                l1 = atv_ijk[Rnh][0];
                l2 = atv_ijk[Rnh][1];
                l3 = atv_ijk[Rnh][2];

                kRn = 2.0*PI*(k2[0]*(double)l1 + k2[1]*(double)l2 + k2[2]*(double)l3);
                kRn -= (dk[0]*Gxyz[ct_AN][0] + dk[1]*Gxyz[ct_AN][1] + dk[2]*Gxyz[ct_AN][2]);

                si1 = sin(kRn);
                co1 = cos(kRn);

                kRn = -2.0*PI*(k1[0]*(double)l1 + k1[1]*(double)l2 + k1[2]*(double)l3);
                kRn -= (dk[0]*Gxyz[ct_AN][0] + dk[1]*Gxyz[ct_AN][1] + dk[2]*Gxyz[ct_AN][2]);

                si2 = sin(kRn);
                co2 = cos(kRn);
                
                for (i1=0; i1<tnoA; i1++){
                    for (j1=0; j1<tnoB; j1++){

                        tmp1r = Wk1[Anum+i1].r*Wk2[Bnum+j1].r + Wk1[Anum+i1].i*Wk2[Bnum+j1].i + Wk1[Anum+fsize+i1].r*Wk2[Bnum+fsize+j1].r + Wk1[Anum+fsize+i1].i*Wk2[Bnum+fsize+j1].i;
                        tmp1i = Wk1[Anum+i1].r*Wk2[Bnum+j1].i - Wk1[Anum+i1].i*Wk2[Bnum+j1].r + Wk1[Anum+fsize+i1].r*Wk2[Bnum+fsize+j1].i - Wk1[Anum+fsize+i1].i*Wk2[Bnum+fsize+j1].r;

                        tmp2r1 = co1*tmp1r - si1*tmp1i;
                        tmp2i1 = co1*tmp1i + si1*tmp1r;

                        tmp1r = Wk1[Bnum+j1].r*Wk2[Anum+i1].r + Wk1[Bnum+j1].i*Wk2[Anum+i1].i + Wk1[Bnum+fsize+j1].r*Wk2[Anum+fsize+i1].r + Wk1[Bnum+fsize+j1].i*Wk2[Anum+fsize+i1].i;
                        tmp1i = Wk1[Bnum+j1].r*Wk2[Anum+i1].i - Wk1[Bnum+j1].i*Wk2[Anum+i1].r + Wk1[Bnum+fsize+j1].r*Wk2[Anum+fsize+i1].i - Wk1[Bnum+fsize+j1].i*Wk2[Anum+fsize+i1].r;

                        tmp2r2 = co2*tmp1r - si2*tmp1i;
                        tmp2i2 = co2*tmp1i + si2*tmp1r;

                        tmp3r = OLPe[bindx][ct_AN][h_AN][i1][j1].r;
                        tmp3i = OLPe[bindx][ct_AN][h_AN][i1][j1].i;

                        sumr += 0.5*(tmp2r1+tmp2r2)*tmp3r - 0.5*(tmp2i1+tmp2i2)*tmp3i;
                        sumi += 0.5*(tmp2r1+tmp2r2)*tmp3i + 0.5*(tmp2i1+tmp2i2)*tmp3r;
                    }
                }
            }   
        }       

        Sop[0] = sumr;
        Sop[1] = sumi;
    }
}