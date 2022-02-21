#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <time.h>
#include <string.h>

#include "mkl.h"
#include "izh.h"
#include "ntk.h"
#include "utils.h"

#define PRINT_PROG
#include "porta_function.h"

extern double _dt;
extern double _R;
extern int print_prog;

// for testing code

int main()
{
    _dt = 0.05;
    _R = 0.05;
    print_prog = 1;

    // init_random_stream(time(NULL));
    init_random_stream(1000);
    double *vm = run_simulation(10, "./data/result");
    free(vm);
    free_rand_stream();
}

