#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <time.h>
#include <string.h>

#include "mpi.h"

#include "mkl.h"
#include "izh.h"
#include "mt64.h"
#include "ntk.h"
#include "utils.h"
#include "porta_function.h"

extern int print_prog;
extern double _dt;
extern double _R;

static void print_elapsed_custom(struct timeval start_t, char *text);

int main(int argc, char **argv)
{
    int world_size, world_rank;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    _dt = 0.05;
    _R = 0.05;
    print_prog = 0;

    double tmax = 10000;
    int max_step = tmax / _dt;

    // init_genrand64(time(NULL) + world_rank*10);
    fprintf(stderr, "Node %d Start\n", world_rank);
    init_random_stream(time(NULL) + world_rank*10);

    // iterate 100 times
    int nitr = 300;
    char fname[100], text[100];
    struct timeval tic;
    

    for (int n=world_rank; n<nitr; n+=world_size)
    {
        gettimeofday(&tic, NULL);

        
        double *vm = run_simulation(10000, NULL);

        sprintf(fname, "./data_mpi/res%d_v.dat", n);
        FILE *fp = fopen(fname, "wb");
        fwrite(vm, sizeof(double), max_step*2, fp);
        fclose(fp);

        free(vm);
        sprintf(text, "itr %02d done,", n);
        print_elapsed_custom(tic, text);
    }

    free_rand_stream();
    MPI_Barrier(MPI_COMM_WORLD);

    MPI_Finalize();

}


static void print_elapsed_custom(struct timeval start_t, char *text)
{
    int sec, msec, usec, x;
    struct timeval end_t;

    gettimeofday(&end_t, NULL);


    sec = end_t.tv_sec - start_t.tv_sec;
    usec = end_t.tv_usec - start_t.tv_usec;
    x = usec / 1e3;
    msec = x;
    usec -= x * 1e3;

    if (usec < 0){
        msec -= 1;
        usec += 1e3;
    }
    if (msec < 0){
        sec -= 1;
        msec += 1e3;
    }

    fprintf(stderr, "%s elapsed time = %ds %dms %dus\n", text, sec, msec, usec);

}
