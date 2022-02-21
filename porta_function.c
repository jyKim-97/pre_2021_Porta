#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <time.h>
#include <string.h>
#include <math.h>

#include "porta_function.h"
#include "mkl.h"
#include "mt64.h"
#include "utils.h"

// #define DEBUG
#define PRINT_NTK
#define PRINT_VAR
#define N_VAR 4

#ifdef DEBUG
    #define print_checkpoint(x) fprintf(stderr, x)
#else
    #define print_checkpoint(x) NULL
#endif

static void write_dat_ind(FILE *fp, int nx, double *x, int id_out[]);
double get_r(int node_post, syn_t *syns, int *cell_types, int pre_type);

/*** PARAMETERS IN THE PAPER ***/

// extern double gRatio;
extern double _R;
extern double default_cell_params[][4];

double cell_ratio[2] = {0.8, 0.2};
double tauAMPA = 5.26;
double tauGABA = 5.6; 
double i_dc = 25; // DC current to pop1

// bck syns
double R[2] = {3000, 2400};
// double R[2] = {2800, 2200};
// double exp(-R[2] = {3000, 2400});
double gP = 0.6;  // for each neuron

// synapse
int num_pre_same = 50;  // type X
int num_pre_other = 20; // from other population's excitatory neuron
// original paper
double gAMPA[2] = {3, 0.8};  // nS
double gGABA[2] = {16, 16.4}; // nS
double gSend[2] = {0.15, 4};  // 1->2, 2->1

int num_cells = 1000;
int num_pop = 500;

/*** ######### END ######### ***/

int print_prog = 0;


// void test()
// {
//     neuron_t cells;
//     syn_t syns;
//     syn_t bck_syns;

//     char fname[200];

//     int *cell_types = (int*) calloc(num_cells, sizeof(int));
//     init_simulation(cell_types, &cells, &syns, &bck_syns);
    
//     int *id_fire = (int*) malloc(sizeof(int) * (num_cells+1));

// }


double *run_simulation(double tmax, char tag[])
{
    neuron_t cells;
    syn_t syns;
    syn_t bck_syns;

    char fname[200];

    int *cell_types = (int*) calloc(num_cells, sizeof(int));
    init_simulation(cell_types, &cells, &syns, &bck_syns);

    int max_step = tmax/_dt;

    progbar_t bar;
    if (print_prog == 1){
        init_progressbar(&bar, max_step);
    }

    struct timeval tic;
    double *vm = (double*) malloc(sizeof(double) * 2* max_step);

    // #ifndef MPI
    // save env
    sprintf(fname, "%s_env.txt", tag);
    FILE *fp = fopen(fname, "w");
    fprintf(fp, "%f-%f\n", tmax, _dt);
    for (int n=0; n<num_cells; n++) { fprintf(fp, "%d,", cell_types[n]); }
    fclose(fp);

    sprintf(fname, "%s_v.txt", tag);
    fp = fopen(fname, "w");

    FILE *fpt;
    sprintf(fname, "%s_spk.txt", tag);
    fpt = fopen(fname, "w");
    // #endif

    #ifdef PRINT_VAR
    int id_out[N_VAR] = {100, 420, 550, 920};
    FILE *ftest[3];
    ftest[0] = fopen("./data/v_out.dat", "wb");
    ftest[1] = fopen("./data/i_out.dat", "wb");
    // ftest[2] = fopen("./data/r_out.dat", "wb");
    #endif

    // // set monitor neurons
    // double **ptr_vm[2];
    // int used[1000] = {0,};
    // for (int i=0; i<2; i++){
    //     ptr_vm[i] = (double**) malloc(sizeof(double*) * 100);
    //     int j=0;
    //     while (j < 100){
    //         int id = genrand64_real2()*400 + 500*i;
    //         if (used[id] == 0){
    //             ptr_vm[i][j] = cells.v + id;
    //             used[id] = 1;
    //             j++;
    //         }
    //     }
    // }

    gettimeofday(&tic, NULL);

    // current 
    double *ic = (double*) calloc(num_cells, sizeof(double));
    for (int n=0; n<num_pop*cell_ratio[0]; n++){
        ic[n] = i_dc;
    }
    int *id_bck_fire = (int*) malloc(sizeof(int) * (num_cells+1));
    for (int n=0; n<max_step; n++){

        memcpy(cells.ic, ic, sizeof(double)*num_cells);

        gen_bck_spike(&bck_syns, id_bck_fire);
        update_bck_porta(&bck_syns, id_bck_fire);
        add_bcksyn_current_porta(&bck_syns, &cells);

        update_syns_porta(&syns, cells.id_fired);
        add_syn_current_porta(&syns, &cells);
 
        update_neurons(&cells, n);
        
        // if (n%100 == 0){
            for (int i=0; i<2; i++){
                vm[2*n+i] = 0;
                for (int j=0; j<400; j++){
                    vm[2*n+i] += cells.v[num_pop*i + j];
                }
                vm[2*n+i] /= 400;
            }
        //     fprintf(fp, "%f,%f,%f,\n", n*_dt, vm[0], vm[1]);
        // }

        // #ifndef MPI
        fprintf(fp, "%f,%f,%f,\n", n*_dt, vm[0], vm[1]);
        write_cell_fire(fpt, cells.id_fired);
        // #endif

        #ifdef PRINT_VAR
        write_dat_ind(ftest[0], N_VAR, cells.v, id_out);
        write_dat_ind(ftest[1], N_VAR, cells.ic, id_out);
        #endif

        if (print_prog == 1){
            progressbar(&bar, n);
        }

    }

    // #ifndef MPI
    fclose(fp);
    fclose(fpt);
    // #endif

    #ifdef PRINT_VAR
    fclose(ftest[0]);
    fclose(ftest[1]);
    #endif

    free(ic);
    free(id_bck_fire);
    free(cell_types);
    free_neurons(&cells);
    free_syns(&syns);
    free_syns(&bck_syns);

    if (print_prog == 1){
        printf("\nDone, ");
        print_elapsed(tic);
    }

    mkl_free_buffers();

    return vm;
}


void gen_celltypes(int *cell_types)
{
    int mod, num_exc=num_pop*cell_ratio[0];
    // 0: pop0 - exc, 1: pop0 - inh, 2: pop1 - exc, 3: pop1 - inh
    for (int n=0; n<num_cells; n++){
        cell_types[n] = n/num_pop*2; // population 0/1
        mod = n%num_pop;
        if (mod > num_exc){ // exc/inh
            cell_types[n]++;
        }
    }
}


void init_simulation(int *cell_types, neuron_t *cells, syn_t *syns, syn_t *bck_syns)
{
    /*** generate cells ***/
    int num_exc=num_pop*cell_ratio[0];
    
    gen_celltypes(cell_types);
    init_cell_vars(cells, num_cells, default_cell_params, cell_types);
    print_checkpoint("Create cell done\n");

    double buf[4];
    for (int n=0; n<num_cells; n++){
        if (cell_types[n]%2 == 0){
            gen_exc_neuron_params(buf);
        } else {
            gen_inh_neuron_params(buf);
        }
        cells->a[n] = buf[0];
        cells->b[n] = buf[1];
        cells->c[n] = buf[2];
        cells->d[n] = buf[3];
    }

    /*** create background synapse ***/
    ntk_t ntk_bck_syns;
    init_bi_ntk(num_cells, num_cells, &ntk_bck_syns);
    memcpy(ntk_bck_syns.node_types, cell_types, sizeof(int)*num_cells);
    for (int n=0; n<num_cells; n++){
        ntk_bck_syns.num_edges[n] = 1;
        ntk_bck_syns.adj_list[n][0] = n;
        ntk_bck_syns.strength[n][0] = gP;
    }

    double bck_veq[4] = {0,};
    double bck_tau[4] = {tauAMPA, tauAMPA, tauAMPA, tauAMPA};
    init_syn_vars(bck_syns, num_cells, BACKGROUND, &ntk_bck_syns, bck_veq, bck_tau, cells->v, cells->ic);
    free_bi_ntk(&ntk_bck_syns);

    for (int n=0; n<num_cells; n++){
        int ctp_exc = n / num_pop;
        // bck_syns->p_fire[n] = R[ctp_exc]/1000*_dt;
        bck_syns->p_fire[n] = 1.-exp(-R[ctp_exc]/1000*_dt);
    }

    /*** create synapse ***/
    ntk_t ntk_syns;
    init_bi_ntk(num_cells, num_cells, &ntk_syns);
    memcpy(ntk_syns.node_types, cell_types, num_cells*sizeof(int));

    // connect within population
    double gbar[4] = {gAMPA[0], gGABA[0], gAMPA[1], gGABA[1]};

    int id_pre0;
    for (int n=0; n<num_cells; n++){
        id_pre0 = (n/num_pop) * num_pop;
        connect_pre(&ntk_syns, n, num_pre_same, id_pre0, id_pre0+num_pop, gbar, cell_types);
    }
    print_checkpoint("connect_within done\n");

    // connect without population
    gbar[0] = gSend[0];
    gbar[2] = gSend[1];
    for (int n=0; n<num_cells; n++){
        id_pre0 = num_pop * (1 - n/num_pop);
        connect_pre(&ntk_syns, n, num_pre_other, id_pre0, id_pre0+num_exc, gbar, cell_types);
    }

    print_checkpoint("connect_without done\n");

    #ifdef PRINT_NTK
    print_adjlist(&ntk_syns, "./data/ntk_adjlist.txt");
    print_strength(&ntk_syns, "./data/ntk_strength.txt");
    #endif

    // gen synapse
    double syn_veq[4] = {0, -65, 0, -65};
    double syn_tau[4] = {tauAMPA, tauGABA, tauAMPA, tauGABA};
    init_syn_vars(syns, num_cells, DELAY, &ntk_syns, syn_veq, syn_tau, cells->v, cells->ic);

    free_bi_ntk(&ntk_syns);
    print_checkpoint("connect synapse done\n");

    #ifdef PRINT_NTK
    print_syn_ntk(syns, "./data/syn_info.txt");
    #endif
}


void connect_pre(ntk_t *ntk, int id_node, int num_pre, int id_pre0, int id_pre1, double gbar[4], int *cell_types)
{
    int n_pick, n_target, n_cum=0;

    n_target = id_pre1 - id_pre0;
    int *picks = (int*) malloc(num_pre * sizeof(int));
    int *used  = (int*) calloc(n_target, sizeof(int));

    while (n_cum < num_pre){
        n_pick = (int) (genrand64_real2() * n_target);

        if (used[n_pick] == 0){
            used[n_pick]   = 1;
            int node_append = n_pick + id_pre0;
            if (node_append != id_node){
                picks[n_cum++] = node_append;
            }
        }
    }
    free(used);


    int *ptr_num_edge;

    for (int n=0; n<num_pre; n++){
        n_pick = picks[n];
        ptr_num_edge = ntk->num_edges + n_pick;

        append_node(id_node, ptr_num_edge, ntk->adj_list+n_pick);
        (*ptr_num_edge)--;

        append_value(gbar[cell_types[n_pick]], ptr_num_edge, ntk->strength+n_pick);
    }
    free(picks);
}


void gen_exc_neuron_params(double buf[4])
{
    // buf: {a, b, c, d}
    double sgm = genrand64_real2();
    double sgm2 = sgm * sgm;

    buf[0] = 0.02;
    buf[1] = 0.20;
    buf[2] = -65 + 15*sgm2;
    buf[3] = 8 - 6*sgm2;
}


void gen_inh_neuron_params(double buf[4])
{
    // buf: {a, b, c, d}
    double sgm = genrand64_real2();

    buf[0] = 0.02 + 0.08 * sgm;
    buf[1] = 0.25 - 0.05 * sgm;
    buf[2] = -65;
    buf[3] = 2;
}


int *get_expand_spk(syn_t *syns, int *id_fire)
{
    int *id_fire_ex = (int*) malloc(sizeof(int) * (syns->num_syns+1));
    int *id_pre = syns->id_pre_neuron;

    int nid = 0;
    int *ptr_id_ex = id_fire_ex;
    int *ptr_id = id_fire;

    while (*ptr_id > -1){
        while (id_pre[nid] - *ptr_id < 0){
            nid++;
        }

        if (id_pre[nid] == *ptr_id){
            *ptr_id_ex++ = nid++;
        } else {
            ptr_id++;
        }
    }
    *ptr_id_ex = -1;

    return id_fire_ex;
}


double *f_dr_syn_porta(double *r, void *arg_syn, void *arg_id_fire)
{
    syn_t *syns = (syn_t*) arg_syn;
    int num_syns = syns->num_syns;

    double *dr = (double*) malloc(sizeof(double) * num_syns);
    memcpy(dr, r, sizeof(double) * num_syns);
    cblas_dscal(num_syns, -_dt, dr, 1);

    int *ptr_id = arg_id_fire;
    while (*ptr_id > -1){
        dr[*ptr_id] += syns->weight[*ptr_id] * _R;
        ptr_id++;
    }

    vdMul(num_syns, syns->inv_tau, dr, dr);

    return dr;
}


void update_syns_porta(syn_t *syns, int *id_fire)
{
    int num_syns = syns->num_syns;

    int *id_fire_ex = get_expand_spk(syns, id_fire);

    double *dr = solve_deq_using_euler(f_dr_syn_porta, num_syns, syns->r, syns, id_fire_ex);
    cblas_daxpy(num_syns, 1, dr, 1, syns->r, 1);
    free(id_fire_ex);
    free(dr);
}


void add_syn_current_porta(syn_t *syns, neuron_t *cells)
{
    int num_syns = syns->num_syns;

    double *isyn = (double*) calloc(num_cells, sizeof(double));
    for (int n=0; n<syns->num_syns; n++){
        int id = syns->id_post_neuron[n];
        isyn[id] += syns->r[n] * (syns->veq[n] - cells->v[id]);
    }
    cblas_daxpy(num_cells, 1, isyn, 1, cells->ic, 1);

    // double *isyn = (double*) malloc(sizeof(double) * num_syns);
    // read_ptr(num_syns, isyn, syns->ptr_vpost);
    // cblas_daxpy(num_syns, -1, syns->veq, 1, isyn, 1); // v - veq
    // vdMul(num_syns, syns->r, isyn, isyn); // r*(v - veq)
    // for (int n=0; n<num_syns; n++){
    //     *(syns->ptr_ipost[n]) -= isyn[n];
    // }
    free(isyn);
}


double *f_dr_bck_porta(double *r, void *arg_syn, void *arg_id_fire)
{
    syn_t *bck_syns = (syn_t*) arg_syn;
    
    double *dr = (double*) malloc(sizeof(double) * num_cells);
    memcpy(dr, r, sizeof(double) * num_cells);
    cblas_dscal(num_cells, -_dt, dr, 1);

    int *ptr_id = arg_id_fire;
    while (*ptr_id > -1){
        dr[*ptr_id] += gP * _R;
        ptr_id++;
    }

    cblas_dscal(num_cells, 1./tauAMPA, dr, 1);
    
    return dr;
}


void update_bck_porta(syn_t *syns, int *id_fire)
{
    double *dr = solve_deq_using_euler(f_dr_bck_porta, num_cells, syns->r, (void*) syns, (void*) id_fire);
    cblas_daxpy(num_cells, 1, dr, 1, syns->r, 1);
    free(dr);
}


void add_bcksyn_current_porta(syn_t *bck_syns, neuron_t *cells)
{
    double *isyn = (double*) malloc(sizeof(double) * num_cells);
    vdMul(num_cells, bck_syns->r, cells->v, isyn);
    cblas_daxpy(num_cells, -1, isyn, 1, cells->ic, 1);
    free(isyn);
}


static void write_dat_ind(FILE *fp, int nx, double *x, int id_out[])
{
    double *x_write = (double*) malloc(nx * sizeof(double));
    int *ptr = id_out;
    for (int n=0; n<nx; n++){
        x_write[n] = x[*ptr++];
    }

    fwrite(x_write, sizeof(double), nx, fp);
    free(x_write);
}


double *solve_deq_using_euler(double* (*f) (double*, void*, void*), int N, double *x, void *arg1, void *arg2)
{
    double *df = f(x, arg1, arg2);
    return df;
}


void write_cell_fire(FILE *fid, int *id_fire)
{
    int *id = id_fire;
    while (*id != -1){
        fprintf(fid, "%d,", *(id++));
    }
    fprintf(fid, "\n");
}


void print_adjlist(ntk_t *ntk, char fname[100])
{
    FILE *fp = fopen(fname, "w");

    for (int n=0; n<ntk->num_pre; n++){
        fprintf(fp, "%d:", n);
        for (int i=0; i<ntk->num_edges[n]; i++){
            fprintf(fp, "%d,", ntk->adj_list[n][i]);
        }
        fprintf(fp, "\n");
    }
    fclose(fp);
}


void print_strength(ntk_t *ntk, char fname[100])
{
    FILE *fp = fopen(fname, "w");

    for (int n=0; n<ntk->num_pre; n++){
        fprintf(fp, "%d:", n);
        for (int i=0; i<ntk->num_edges[n]; i++){
            fprintf(fp, "%f,", ntk->strength[n][i]);
        }
        fprintf(fp, "\n");
    }
    fclose(fp);
}


void print_syn_ntk(syn_t *syns, char fname[100])
{
    FILE *fp = fopen(fname, "w");
    for (int n=0; n<syns->num_syns; n++){
        fprintf(fp, "%d-%d-%f\n", syns->id_pre_neuron[n], syns->id_post_neuron[n], syns->weight[n]);
    }
    fclose(fp);
}

