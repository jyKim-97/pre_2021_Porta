#!/bin/bash

if [ "$1" == "mpi" ]; then
    CFLAGS="mpicc"
    TARGET="main_mpi.c"
    OUT="main_mpi.out"
    ENV="-DMPI"
elif [ "$1" == "test" ]; then
    CFLAGS="gcc"
    TARGET="test.c"
    OUT="test.out"
    ENV=""
else
    CFLAGS="gcc"
    TARGET="main.c"
    OUT="main.out"
    ENV=""
fi

echo "$TARGET -> $OUT"

OPT="-g -O3 -Wall -std=c11"
INC="-I${MKLROOT}/include -I${IZHROOT}/include"
LIB="-L${MKLROOT}/lib/intel64 -L${IZHROOT}/lib"
LDFLAGS="-lizh_neuralnet -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lm"

rm -f $OUT
$CFLAGS $OPT -o tmp.o -c $TARGET $INC $LIB $LDFLAGS
$CFLAGS $OPT $ENV -o ./porta_function.o -c ./porta_function.c $INC $LIB $LDFLAGS
$CFLAGS $OPT -o $OUT ./porta_function.o ./tmp.o $INC $LIB $LDFLAGS
