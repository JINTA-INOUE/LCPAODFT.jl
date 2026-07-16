/usr/local/bin/gcc-15 -O3 -shared -fPIC -o Generate_MLWF.so \
    -llapack -lopenblas \
    -I/usr/local/opt/openblas/include -L/usr/local/opt/openblas/lib \
    sgn.c EigenState_k.c Projection_Amatrix.c Calc_Mmnkb_zero.c Write_MLWF.c \
    Generate_MLWF.c -ffast-math -lm -fopenmp