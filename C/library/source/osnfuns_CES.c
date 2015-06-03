#include "osnfuns_CES.h"

/*======================================================*
 * FUNCTION: Calc_Dimensionality						*
 * 		Calculates the number of one-step neighbors for *
 * 		a given Seq										*
 * =====================================================*/

int Calc_Dimensionality(int *codon_index_vec, int aa_count)
{
	int i, Dim;

	Dim = 0;
	for (i = 1; i < aa_count; i++) {
		Dim += Codon[codon_index_vec[i]].num_one_step_synonym;
	}

	return Dim;

}

/*==================================================================================*
 * FUNCTION: Generate_D_Array							                            *
 * 	     D_array is a 2D matrix of size 2 x Dimensionality of sequence              *
 *              Dimensionality = # one step syonymous neighbors                     *
 *           D_array[i][0] = j, where j is the position of an amino acid in the ORF *
 *           D_array[i][1] = k, where k is a codon_index value for an alternative   *
 *                            synonymous codon for the wt sequence                  *
 *                                                                                  *
 * =================================================================================*/

int
Generate_D_Array(int *codon_index_vec,
		 int *D_array[][2], int aa_count)
{
	int i, j, jmax;
	int D_index;
	int codon_index;

	//  imax = (int)(aa_count);

	D_index = 0;
	for (i = 1; i < aa_count; i++) {
		codon_index = codon_index_vec[i];
		jmax = Codon[codon_index].num_one_step_synonym;

		//go through each synonym
		for (j = 0; j < jmax; j++) {
			*(D_array[D_index][0]) = (int)(i);
			*(D_array[D_index++][1]) =
			    (Codon[codon_index].one_step_synonym_index[j]);
		}

	}
	return 0;
}

/*==========================================================*
 * FUNCTION: Calc_Delta_Eta_NSE_First_and_Second_Term_Vecs	*
 * 		Used in conjunction with Calc_Delta_Eta_NSE to		*
 * 		calculate delta eta values for all one step 		*
 * 		neighbors											*
 * =========================================================*/

void
Calc_Delta_Eta_NSE_First_and_Second_Term_Vecs(double *sigma_ratio_vec,
					      double *b_over_c_vec, double
					      *delta_eta_first_term_vec, double
					      *delta_eta_second_term_vec,
					      int aa_count)
{
	int i;
	double curr_cost;
	double delta_eta_first_term;	//, delta_eta_second_term;

	curr_cost = (double)(A1);

	*(delta_eta_first_term_vec++) = 0.0;	//zero b/c there is no chance of a nonsense error at the first codon.
	*(delta_eta_second_term_vec++) = 0.0;	//zero b/c there is no chance of a nonsense error at the first codon.

	b_over_c_vec++;
	sigma_ratio_vec++;

	delta_eta_first_term = 0.;

	for (i = 1; i < aa_count; i++) {
		//term is a summation from j=0 to i-1, so we don't update the first term until the end
		*(delta_eta_first_term_vec++) = delta_eta_first_term;

		curr_cost += A2;

		*(delta_eta_second_term_vec) =
		    curr_cost * (*(sigma_ratio_vec++));

		delta_eta_first_term +=
		    (*(b_over_c_vec++)) * (*(delta_eta_second_term_vec++));

		//(A1 + A2 i)*b_over_c_vec[i]*sigma_ratio_vec[i];
	}

}

/*===========================================================*
 * FUNCTION: Calc_Delta_Eta_NSE_Vec						     *
 * 		Used in conjunction with Calc_Delta_Eta_NSE_First... *
 * 		to calculate delta eta values for all one step 	     *
 * 		neighbors										     *
 * ==========================================================*/

void
Calc_Delta_Eta_NSE_Vec(int *codon_index_vec,
		       double *b_over_c_vec, double *first_term_vec,
		       double *second_term_vec, double *delta_eta_vec,
		       int aa_count, int D_dim)
{
	int i, j, jmax, D_index;
	double curr_cost;
	int codon_index;
	double first_term, second_term;
	double z;		//z = \frac{\?igma_{n,i}}{\sigma_{n,j}} =  \left(\frac{c_{k,i}}{c_{k,i}+b} \frac{c_{k,j}+b}{c_{k,j}}\right)

	D_index = 0;

	curr_cost = (double)(A1);
	codon_index_vec++;
	first_term_vec++;
	second_term_vec++;

	//go through each codon
	for (i = 1; i < aa_count; i++) {
		codon_index = *(codon_index_vec++);
		jmax = Codon[codon_index].num_one_step_synonym;

		curr_cost += A2;

		first_term = *(first_term_vec++);
		second_term = *(second_term_vec++);
		//go through each synonym
		for (j = 0; j < jmax; j++) {
			z = Codon[codon_index].elong_pr_ratio[j];

			delta_eta_vec[D_index++] =
			    first_term * (1 - z) +
			    second_term * Codon[codon_index].delta_b_over_c[j];
		}

	}
	if (D_index != D_dim) {	//not sure why I need to add one before comparing to D_dim
		if(R_or_C){
			printf("ERROR in Calc_Delta_Eta_NSE_Vec: %d out of %d indices cycled through. Exiting",
				 D_index, D_dim);
		}else {
			Rprintf("ERROR in Calc_Delta_Eta_NSE_Vec: %d out of %d indices cycled through. Exiting",
				 D_index, D_dim);
	    }
		exit(1);
	};
}

/*======================================================*
 * FUNCTION: Calc_Pi_Vec								*
 * 		Calculates pi(i->j) value according to 			*
 * 		Sella and Hirsh and stores in pi_vec. 			*
 * =====================================================*/

void
Calc_Pi_Vec(double *delta_eta_vec, double *pi_vec, double *ptr_pi_total,
	    double phi, double qPhi, double two_qNePhi, int aa_count, int D_dim)
{
	int i;
	double pi, delta_eta;
	double pi_total;
	double invNe;

	//assuming diploid Wright-Fischer Process.  If haploid need to multiply qPhi by 2
	invNe = 1. / Ne;

	pi_total = 0;
	for (i = 0; i < D_dim; i++) {
		delta_eta = delta_eta_vec[i];

		if (delta_eta == 0 || phi == 0) {	//was producing strange behavior when phi=0.  First codon was always changing.
			//pure drift process
			pi = pi_vec[i] = invNe;
		} else {
			//numerator = -expm1(minuxtwoqphi*delta_eta);//expm1(x) = exp(x)-1, -expm1 = (1-exp(x))
			//denominator = -expm1(minustwoqNephi*delta_eta);
			pi = pi_vec[i] = expm1(-qPhi * delta_eta) / expm1(-two_qNePhi * delta_eta);	//numerator/denominator;
		}
		pi_total += pi;
	}
	*(ptr_pi_total) = pi_total;
}

/*======================================================*
 * FUNCTION: Calc_Replacemet_Pr_Vec						*
 * 		Weights pi by mu to get replacment pr			*
 * =====================================================*/

void
Calc_Replacement_Pr_Vec(double *mutation_vec, double *pi_vec, double *pr_vec,
			double *ptr_pr_total, int D_dim)
{
	int i;
	double pr;
	double pr_total;

	pr_total = 0;
	for (i = 0; i < D_dim; i++) {
		pr = pr_vec[i] = mutation_vec[i] * pi_vec[i];
		pr_total += pr;
	}

	*(ptr_pr_total) = pr_total;
}
