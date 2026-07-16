#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "LCPAO2Wannier.h"

void Set_Gmnk2rimnk(int Nkpt, int WANNUM, dcomplex *rimnk, dcomplex ***Gmnk)
{
    int bindx, k, j, i;

    bindx=0;
    for(k=0;k<Nkpt;k++){
        for(j=0;j<WANNUM;j++){
            for(i=0;i<WANNUM;i++){
                rimnk[bindx].r=Gmnk[k][i][j].r;
                rimnk[bindx].i=Gmnk[k][i][j].i;
                bindx++;
            }
        }
    }
}