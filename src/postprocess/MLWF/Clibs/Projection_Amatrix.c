#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "Generate_MLWF.h"


#pragma optimization_level 1
void Projection_Amatrix(
    int MLWF_Num_Kinds_Projectors, int **MLWF_NumL_Pro,
    int *MLWF_Num_Pro, int **MLWF_Select_Matrix, double ***MLWF_Projector_Hybridize_Matrix, double ****MLWF_RotMat_for_Real_Func,
    int *FNAN_WP, int **natn_WP, int **ncn_WP, int **atv_ijk, int *Total_NumOrbs, 
    double ****OLP_WP,
    dcomplex ****Amnk, double **kg, int spinsize, 
    int fsize, int SpinP_switch, 
    int Nkpt, int band_num,  int wan_num, 
    dcomplex ****Wkall, int *MP, int ***Nk)
{
    /* calculate the A matrix  */   
    int ct_AN, h_AN, Gh_AN, i,j,k, TNO1, TNO2, Nspink;
    int spin, Rn, disentangle, BAND;
    int nindx, windx, proj_kind, tot_loc_basis;
    int mu1, L, tnoB, Anum, Bnum, Rnh, l1, l2, l3, i1,p;
    double co, si, kRn, sumr, sumi, tmpr, tmpi;
    double sr, pxr, pxi, pyr, pyi, pzr, pzi;
    dcomplex **tmpAmn; /* temporary Amn matrix*/ 
    dcomplex *tmpResult;



    /* count totally how many local orbitals are included */
    tot_loc_basis=0;
    for(i=0;i<MLWF_Num_Kinds_Projectors;i++){
        for(L=0;L<=3;L++){
            tot_loc_basis+=MLWF_NumL_Pro[i][L]*(2*L+1);
        }
    } 

    if(tot_loc_basis >=(2*3+1)){
        if(tot_loc_basis>=wan_num){
            tmpResult=(dcomplex*)malloc(sizeof(dcomplex)*tot_loc_basis);
        }
        else{
            tmpResult=(dcomplex*)malloc(sizeof(dcomplex)*wan_num);
        }
    }else{
        if(wan_num<2*3+1){
            tmpResult=(dcomplex*)malloc(sizeof(dcomplex)*(2*3+1));
        }else{
            tmpResult=(dcomplex*)malloc(sizeof(dcomplex)*wan_num);
        }
    }

    tmpAmn=(dcomplex**)malloc(sizeof(dcomplex*)*band_num);
    for(i=0;i<band_num;i++){
        tmpAmn[i]=(dcomplex*)malloc(sizeof(dcomplex)*tot_loc_basis);
    }

    BAND=band_num;
    if(band_num>wan_num){
        disentangle=1;
    }else{
        disentangle=0;
    }
    

    for(spin=0;spin<spinsize;spin++){
        for(k=0;k<Nkpt;k++){
            if(disentangle){
                band_num=Nk[spin][k][1]-Nk[spin][k][0];
            }else{
                band_num=wan_num;
            }

            Nspink = Nk[spin][k][0]+1;

            for(mu1=0;mu1<band_num;mu1++){
                windx=0;
                for(proj_kind=0;proj_kind<MLWF_Num_Kinds_Projectors; proj_kind++){
                    Anum = 0;
                    for (L=0; L<=3; L++){
                        Anum += (2*L+1)*MLWF_NumL_Pro[proj_kind][L];
                    }

                    for(nindx=0;nindx<Anum;nindx++){

                        sumr=0.0;
                        sumi=0.0; 

                        for(h_AN=0;h_AN<FNAN_WP[proj_kind];h_AN++){
                            Rnh=ncn_WP[proj_kind][h_AN];
                            Gh_AN=natn_WP[proj_kind][h_AN]-1;
                            tnoB=Total_NumOrbs[Gh_AN];
                            Bnum=MP[Gh_AN]; 
                            l1=atv_ijk[Rnh][0];
                            l2=atv_ijk[Rnh][1];
                            l3=atv_ijk[Rnh][2];
                            kRn=-2.0*PI*(kg[k][0]*(double)l1+kg[k][1]*(double)l2+kg[k][2]*(double)l3);
                            co=cos(kRn);
                            si=sin(kRn);

                            for(i1=0;i1<tnoB;i1++){
                                if(SpinP_switch!=3){
                                    tmpr = co*Wkall[k][spin][Nspink+mu1][Bnum+i1].r+si*Wkall[k][spin][Nspink+mu1][Bnum+i1].i;
                                    tmpi =-Wkall[k][spin][Nspink+mu1][Bnum+i1].i*co+si*Wkall[k][spin][Nspink+mu1][Bnum+i1].r;
                                }else{
                                    if(windx<tot_loc_basis/2){
                                        tmpr= co*Wkall[k][0][Nspink+mu1][Bnum+i1].r+si*Wkall[k][0][Nspink+mu1][Bnum+i1].i;
                                        tmpi=-Wkall[k][0][Nspink+mu1][Bnum+i1].i*co+si*Wkall[k][0][Nspink+mu1][Bnum+i1].r;
          
                                    }else{
                                        tmpr= co*Wkall[k][0][Nspink+mu1][Bnum+fsize+i1].r+si*Wkall[k][0][Nspink+mu1][Bnum+fsize+i1].i;
                                        tmpi=-Wkall[k][0][Nspink+mu1][Bnum+fsize+i1].i*co+si*Wkall[k][0][Nspink+mu1][Bnum+fsize+i1].r;
                                    }
                                }
                                tmpr=tmpr*OLP_WP[proj_kind][h_AN][nindx][i1];
                                tmpi=tmpi*OLP_WP[proj_kind][h_AN][nindx][i1];

                                sumr=sumr+tmpr;
                                sumi=sumi+tmpi;
                            }
                        }

                        tmpAmn[mu1][windx].r=sumr/sqrt(fabs((double)FNAN_WP[proj_kind])); 
                        tmpAmn[mu1][windx].i=sumi/sqrt(fabs((double)FNAN_WP[proj_kind]));

                        windx++;
                    }
                }
            } 

            for(mu1=0;mu1<band_num;mu1++){
                windx=0;
                for(proj_kind=0;proj_kind<MLWF_Num_Kinds_Projectors; proj_kind++){
                    for (L=0; L<=3; L++){
                        if(MLWF_NumL_Pro[proj_kind][L]!=0 && L!=0){
                            for(i=0;i<2*L+1;i++){
                                sumr=0.0;sumi=0.0;
                                for(j=0;j<2*L+1;j++){
                                    sumr=sumr+MLWF_RotMat_for_Real_Func[proj_kind][L][i][j]*tmpAmn[mu1][windx+j].r;
                                    sumi=sumi+MLWF_RotMat_for_Real_Func[proj_kind][L][i][j]*tmpAmn[mu1][windx+j].i;
                                }
                                tmpResult[i].r=sumr; 
                                tmpResult[i].i=sumi;
                            }
                            for(i=0;i<2*L+1;i++){
                                tmpAmn[mu1][windx+i].r=tmpResult[i].r; 
                                tmpAmn[mu1][windx+i].i=tmpResult[i].i;
                            }
                            windx+=2*L+1;
                        }else if(MLWF_NumL_Pro[proj_kind][L]!=0 && L==0){
                            windx++;
                        }
                    }
                }
            } 

            for(mu1=0;mu1<band_num;mu1++){
                windx=0; 
                nindx=0;
                for(proj_kind=0;proj_kind<MLWF_Num_Kinds_Projectors; proj_kind++){
                    for(i=0;i<MLWF_Num_Pro[proj_kind];i++){
                        Amnk[spin][k][mu1][nindx]=tmpAmn[mu1][MLWF_Select_Matrix[proj_kind][i]+windx]; 
                        nindx++;
                    }
                    for (L=0; L<=3; L++){
                        windx += (2*L+1)*MLWF_NumL_Pro[proj_kind][L];
                    }
                }
            }

            for(mu1=0;mu1<band_num;mu1++){
                nindx=0;
                for(proj_kind=0;proj_kind<MLWF_Num_Kinds_Projectors; proj_kind++){
                    for(i=0;i<MLWF_Num_Pro[proj_kind];i++){
                        sumr=0.0;
                        sumi=0.0;
                        for(j=0;j<MLWF_Num_Pro[proj_kind];j++){ 
                            sumr=sumr+MLWF_Projector_Hybridize_Matrix[proj_kind][i][j]*Amnk[spin][k][mu1][nindx+j].r;
                            sumi=sumi+MLWF_Projector_Hybridize_Matrix[proj_kind][i][j]*Amnk[spin][k][mu1][nindx+j].i;
                        }
                        tmpResult[i].r=sumr;
                        tmpResult[i].i=sumi; 
                    }
                    for(i=0;i<MLWF_Num_Pro[proj_kind];i++){
                        Amnk[spin][k][mu1][nindx+i].r=tmpResult[i].r;
                        Amnk[spin][k][mu1][nindx+i].i=tmpResult[i].i;
                    }

                    nindx+=MLWF_Num_Pro[proj_kind];
                }           
            }
        }
    }


    for(i=0;i<band_num;i++){
        free(tmpAmn[i]); 
    }
    free(tmpAmn);
    free(tmpResult);
}