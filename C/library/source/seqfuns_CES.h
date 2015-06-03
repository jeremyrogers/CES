#ifndef SEQFUNS_CES_HEADER_
#define SEQFUNS_CES_HEADER_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>		//gettimeofday()
#include <gsl/gsl_rng.h>	//uniform rng
#include <gsl/gsl_randist.h>	//rng from distns
#include <math.h>
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

/*=====================*
 * Function Prototypes *
 *=====================*/

void Convert_Codon_Index_to_AA_Index(int *codon_index_vec,
				          int *aa_index_vec, int aa_count);
void Generate_Random_Codon_Seq_and_Index (char codon_seq[][4],
				              int *codon_index,
				              int *aa_index, int aa_count);
void Generate_Random_Codon_Index(int *codon_index,
				    int *aa_index, int aa_count);

void Calc_B_over_C_Vec(int *codon_index_vec, double *b_over_c_vec,
		  int aa_count);

void Convert_Codon_Index_to_Sigma_Vec(int *index_vec,
	                                 double *sigma_vec, int aa_count);
double Calc_Xi_NSE (double *sigma_vec, int aa_count);
double Calc_Eta_NSE (struct seq_struct * seq);

double Calc_Seq_Mu(int * codon_counts_vec, int aa_count);

void Convert_Sigma_Vec_to_Sigma_Ratio_Vec (double *sigma_vec,
				      double *sigma_ratio_vec, int aa_count);
void Recalc_Sigma_Vec (double *sigma_vec, int mut_codon_pos, double factor,
		  int aa_count);


#endif
