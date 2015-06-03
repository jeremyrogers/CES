#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>		//gettimeofday()
#include <R.h>

#include "CES_struct.h"
#include "globvars.h"


/*======================*
 * External definitions *
 *======================*/

extern struct amino_acid AA[22];
extern struct codon_struct Codon[64];
/* We are dealing with an array of sequences in *
 * data simulation, but only one sequence in    *
 * sequence evolution                           */
extern struct seq_struct Seq[6000];
extern struct seq_struct Sequence;

extern double Calc_Eta_NSE(double sigma_n, double xi, int aa_count);

/*=====================*
 * Function Prototypes *
 *=====================*/

int Calc_Dimensionality(int *codon_index_vec, int aa_count);
int Generate_D_Array(int *codon_index_vec,
		         int *D_array[][2], int aa_count);
void Calc_Delta_Eta_NSE_First_and_Second_Term_Vecs(double *sigma_ratio_vec,
		         double *b_over_c_vec, double *delta_eta_first_term_vec, 
		         double *delta_eta_second_term_vec, int aa_count);
void Calc_Delta_Eta_NSE_Vec(int *codon_index_vec,
		         double *b_over_c_vec, double *first_term_vec,
		         double *second_term_vec, double *delta_eta_vec,
		         int aa_count, int D_dim);
void Calc_Pi_Vec(double *delta_eta_vec, double *pi_vec, double *ptr_pi_total,
	             double phi, double qPhi, double two_qNePhi, int aa_count, 
	             int D_dim);
void Calc_Replacement_Pr_Vec(double *mutation_vec, double *pi_vec, double *pr_vec,
			     double *ptr_pr_total, int D_dim);



