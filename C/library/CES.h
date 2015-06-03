#ifndef CES_Header_
#define CES_Header_

#include "CES_struct.h"

/*======================*
 * External Definitions *
 *======================*/

/* FCNS in iofuns_CES.c */

extern int Read_tRNA_File(char *filename);
extern int Process_AA_Information();
extern void Generate_Codon_Structures();
extern int wrong();

/* FCNS in seqfuns_CES.c */

extern void Convert_Codon_Index_to_AA_Index(int *codon_index_vec,
				          int *aa_index_vec, int aa_count);			          
extern void Generate_Random_Codon_Seq_and_Index (char codon_seq[][4],
				              int *codon_index,
				              int *aa_index, int aa_count);
extern void Generate_Random_Codon_Index(int *codon_index,
				    int *aa_index, int aa_count);
extern void Convert_Codon_Index_to_Sigma_Vec(int *index_vec,
	                                 double *sigma_vec, int aa_count);
extern void Calc_B_over_C_Vec(int *codon_index_vec, 
                              double *b_over_c_vec, int aa_count);
extern void Convert_Sigma_Vec_to_Sigma_Ratio_Vec (double *sigma_vec,
				      double *sigma_ratio_vec, int aa_count);
extern void Recalc_Sigma_Vec (double *sigma_vec, int mut_codon_pos, double factor,
		  int aa_count);
extern double Calc_Xi_NSE (double *sigma_vec, int aa_count);
extern double Calc_Eta_NSE (struct seq_struct * seq);
extern double Calc_Seq_Mu(int * codon_counts_vec, int aa_count);

/* FCNS in osnfuns_CES.c */

extern int Calc_Dimensionality(int *codon_index_vec, int aa_count);
extern int Generate_D_Array(int *codon_index_vec,
		         int *D_array[][2], int aa_count);
extern void Calc_Delta_Eta_NSE_First_and_Second_Term_Vecs(double *sigma_ratio_vec,
		         double *b_over_c_vec, double *delta_eta_first_term_vec, 
		         double *delta_eta_second_term_vec, int aa_count);
extern void Calc_Delta_Eta_NSE_Vec(int *codon_index_vec,
		         double *b_over_c_vec, double *first_term_vec,
		         double *second_term_vec, double *delta_eta_vec,
		         int aa_count, int D_dim);
extern void Calc_Pi_Vec(double *delta_eta_vec, double *pi_vec, double *ptr_pi_total,
	             double phi, double qPhi, double two_qNePhi, int aa_count, 
	             int D_dim);
extern void Calc_Replacement_Pr_Vec(double *mutation_vec, double *pi_vec, double *pr_vec,
			     double *ptr_pr_total, int D_dim);

#endif
