#include "f77func.h"

#define PI              3.1415926535897932384626
#define eV2Hartree 27.2113845
#define LAPACK_ABSTOL     6.0e-15    /* absolute error tolerance for lapack routines */
#define smallvalue   1.0e-6
#define BohrR_Wannier    0.529177249
#define BUFFSIZE    2048

#ifndef ___INTEGER_definition___
typedef int INTEGER; /* for fortran integer */
#define ___INTEGER_definition___ 
#endif

#ifndef ___logical_definition___
typedef long int logical;
#define ___logical_definition___ 
#endif

#ifndef ___dcomplex_definition___
typedef struct { double r,i; } dcomplex;
#define ___dcomplex_definition___ 
#endif

double sgn(double nu);
void Write_win(
    char *filename, char **Atoms_Symbol, int Natom, double **Gxyz, double **tv, double **rtv, int BANDNUM, int Wannier_Func_Num, 
    double Wannier_Outer_Window_Bottom, double Wannier_Outer_Window_Top, double Wannier_Inner_Window_Bottom, double Wannier_Inner_Window_Top,
    int knum_i, int knum_j, int knum_k, double **kg);
void Write_eig(char *filename, int BANDNUM, int Nkpt, int spinsize, int fsize3, int ***Nk, double ***EigenValall, double ChemP);
void Write_Amnk(char *filename, int BANDNUM, int Nkpt, int WANNUM, int spinsize, dcomplex ****Amnk);
void Write_Mmnkb(
    char *filename, 
    int BANDNUM, int Nkpt, int tot_bvector, int spinsize,
    double **kg, double **frac_bv, int **kplusb, dcomplex *****Mmnkb_zero);
void Write_Uk(
    char *filename, 
    int BANDNUM, int Nkpt, int WANNUM, int spinsize,
    double **kg, dcomplex ****Uk);
void EigenState_k(
    double k1, double k2, double k3, int *MP, int spinsize, int SpinP_switch, 
    int Natom, int *FNAN, int **natn, int **ncn, int **atv_ijk, int *Total_NumOrbs,
    int fsize, int fsize2, int fsize3,  dcomplex ***Wk1, double **EigenVal1, 
    double ****OLP, double *****Hks, double *****iHks);
void Calc_Mmnkb_zero(
    double k1[3], double k2[3], double dk[3], int bindx,
    int SpinP_switch, int Natom, double **Gxyz, int *FNAN, int **natn, int **ncn, 
    int **atv_ijk, int *Total_NumOrbs, int *MP, int fsize, double Sop[2], 
    dcomplex *****OLPe, dcomplex *Wk1, dcomplex *Wk2);
void Projection_Amatrix(
    int MLWF_Num_Kinds_Projectors, int **MLWF_NumL_Pro,
    int *MLWF_Num_Pro, int **MLWF_Select_Matrix, double ***MLWF_Projector_Hybridize_Matrix, double ****MLWF_RotMat_for_Real_Func,
    int *FNAN_WP, int **natn_WP, int **ncn_WP, int **atv_ijk, int *Total_NumOrbs, 
    double ****OLP_WP,
    dcomplex ****Amnk, double **kg, int spinsize, 
    int fsize, int SpinP_switch, 
    int Nkpt, int band_num,  int wan_num, 
    dcomplex ****Wkall, int *MP, int ***Nk);
