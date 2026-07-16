#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "Generate_MLWF.h"


double sgn(double nu)
{
    if (nu<0.0){
        return -1.0;
    }else{
        return 1.0;
    }
}