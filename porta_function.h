#include "izh.h"
#include "ntk.h"

#ifndef _porta
#define _porta

double *run_simulation(double tmax, char tag[]);
void gen_celltypes(int *cell_types);
void init_simulation(int *cell_types, neuron_t *cells, syn_t *syns, syn_t *bck_syns);
void connect_pre(ntk_t *ntk, int id_node, int num_pre, int id_pre0, int id_pre1, double gbar[4], int *cell_types);
void gen_exc_neuron_params(double buf[4]);
void gen_inh_neuron_params(double buf[4]);

// update codes
int *get_expand_spk(syn_t *syns, int *id_fire);
double *f_dr_syn_porta(double *r, void *arg_syn, void *arg_syn_act);
double *f_dr_bck_porta(double *r, void *arg_syn, void *id_fire);
void update_cells_porta(int nstep, neuron_t *cells);
void update_bck_porta(syn_t *syns, int *id_fire);
void update_syns_porta(syn_t *syns, int *id_fire);
void add_syn_current_porta(syn_t *syns, neuron_t *cells);
void add_bcksyn_current_porta(syn_t *bck_syns, neuron_t *cells);
double *alloc_spk(syn_t *syns, int *id_fire);
double *solve_deq_using_euler(double* (*f) (double*, void*, void*), int N, double *x, void *arg1, void *arg2);
void write_cell_fire(FILE *fid, int *id_fire);
void print_adjlist(ntk_t *ntk, char fname[100]);
void print_strength(ntk_t *ntk, char fname[100]);
void print_syn_ntk(syn_t *syns, char fname[100]);

#endif