/*===========================*
 *  Declaring Header Files   *
 *===========================*/

#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>       //for mkdir
#include <sys/types.h>      //for types w/in mkdir
#include <math.h>
#include <time.h>
#include <sys/time.h>       //gettimeofday()
#include <string.h>
#include <omp.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_min.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_rng.h>    //uniform rng
#include <gsl/gsl_randist.h>    //rng from distns
#include <R.h>

/*==================================*
 * Include CES library header files *
 *==================================*/

#include "CES.h"

/*=======================*
 * Turn off error codes  *
 *=======================*/

int gsl_errors_off = 0;
int print_debug = 0; //1 gives some debug info, 2 gives more

/*========================================================*
 *                DECLARING GLOBAL VARIABLES              *
 *========================================================*/

double A1;          //Cost of initiation of translation in units of ATP. Default val set in command line to  A1 = 2;
double A2;          //Cost each elongation step in ATPs, default set in command line to A2=4
double B;           // Nonsense error rate. units of 1/sec. Default value set in commandline function to   B = 0.00515
double pi = M_PI;       //
double mu = MU;         //per nt mutation rate.
double min_z_max = 1e-10;   //minimum value to start integration routines with;
double max_z_max = 1e100;   //set an upper bound for z
double z_max_factor = 10;   //amount to multiply z_max by when searching for value
double relerr = 1e-6;       //relative error for integration routines
double abserr = 1e-20;      //absolute error for integration routines.  aka tolerance
double at_bias = 0.5;       //parameter for adjusting nt composition due to biased mutation and/or gene conversion.
                            //According to table 6.1 in Lynch (2007) on p. 125, the observed AT bias is 0.62
double mu_bias = 1;         //parameter calculated based on at_bias.  mu_bias = at_bias/(1-at_bias)
double gamma_ratio = 1;     //ratio of transition mutation rate to transversion mutation rate 
                            //This is \alpha/\beta in the Tamura and Nei (1993) where we assume \alpha_1 = \alpha_2
                            //Calculations using data from Lynch give value of approximately 1.22 for S.c. 
double mutation_matrix[4][4];   //mutation matrix order is ATCG.  See mutation.matrices.nb for details.

double global_max_time = MAX_TIME;  //Max number of evolutionary steps to simulate, if set to a negative value it will set the max time so that, under neutrality, there are on average X substitutions/nt.
double Ne = 1;//NE;
double Q = 1;//4.19e-7;     // Fitness scaling coefficient

int ignore_aa;              //decrease true aa_count by this amount.
int random_start = 0;       //indicate whether or not to start with a random sequence (1) or the seq read in (0);
int n_aa = 0,               //Number of aa in a sequence
    n_seq,                  //set when reading fasta file
    pout,                   //print out configuration flag.  Stores -W arguement passed on command line
    pconf, gconf,
    //  calcamean=0,        //indicates whether to calculate arithmetic mean.  Will be reset based on arguments passed
    //  calcgmean=0,
    //  calcvar=0,
    //  calcmode = 1,       //indicates which method to use to estimate the posterior mode. 1:gamma based 2: numerical search
    runs, analytic, max_aa, MLE;
int print_evol_fasta = 0;       //indicate whether or not to print seq evolution to fasta file
                                //*cannot* be set via command line
int print_delta_eta_vals = 0;   //indicate whether or not to print list of delta eta vals for each step
                                //can be set to 1 via command line
int print_eta_trace = 0;        //indicate whether or not to print the eta values and residence times for the wild type genotypes
                                //can be set to 1 via command line

int benchmark = 0;              //Uses predetermined values to generate random sequence and calculate mut_pr and time_step
                                //instead of random generated 
int codon_counts = 0;
int random_walk = 0;

/*=========================================*
 *     Code Specific Global Variables      *
 *  These variables are different in the   *
 *  data simulation and sequence evolution *
 *  code.                                  *
 *=========================================*/

int simulation_steps;
int burn_in_steps;
int total_evol_steps;

int non_random_time_step = 1;   //If this is equal to 0, time step will be pulled from a random exponential dist, if this 
                                //is equal to 1, time step will be calculated based on leaving rate **DEFAULT 1**
int mutation = 1;

int R_or_C = 0; //This flag is used to determine whether the CES library functions should use printf or Rprintf
                //R_or_C = 0 -> Calling from R, use Rprintf

/*========================================================*
 *                 DEFAULT INPUT FILES                    *
 *========================================================*/

char def_tRNA[150] = "tRNA_link";
char fasta_file[150] = "fasta_link";
char out_prefix[150], trna[150], mutation_file[150] = "NONE";
char out_folder[150];       //derived from out_prefix
char command_line[1000];    //save command line for printing
char exec_time[360];        //save time command executed at.

time_t start;
struct tm *timeinfo;
struct timeval time_val, comp_start, comp_stop, print_stop, time_diff;

/*================*
 * Define Structs *
 *================*/

struct amino_acid AA[22];
struct codon_struct Codon[64];
struct seq_struct Sequence;

/*======================================================*
 *              Function Declarations                   *
 *======================================================*/
 
 void Process_R_input(int *CODON_INDEX, double *ELONGATION_RATES, 
             double *MUTATION_RATES, int *AA_COUNT, double *PHI, 
             double *POP, double *A_1, double *A_2,
             double *AT_BIAS, double *BEE, 
             int *IGNORE, double *GAMMA,
             double *SCALE_FACTOR, char ** AA_VEC,
             char ** CODON_VEC, int * N_CODON);
             
 void Rprint_codon_sequence(int *codon_index);

 void Codon_Counts(int * codon_index_vec ,int * codon_cts_vec ,int aa_count);
 
 void bin_eta(double eta, double * bins, double * bin_lims, double value, int num_bins);
 
 void intervals(double from, double to, int length, double * int_lims);
 
 void flip_codons(struct seq_struct *currSeq, struct seq_struct *propSeq, double *randNums, int num_codons);
 
 
/*======================================================*
 * FUNCTION: Evolve_Sequence                            *
 *      This is where the main calculations occur       *
 * =====================================================*/

int
Evolve_Sequence(struct seq_struct *seq)
{
	
	int aa_count = seq->aa_count;
	int aa_position,wt_index,mut_index;
	
    char wt_codon[4];
    char mut_codon[4];

    int D_array[MAX_DDIM][2];    //array with codon position and codon index of one step neighbors
    int mut_codon_index;
    int wt_codon_index;
    int mut_codon_pos;   //position numbers for the codon and nt that's mutated

    int D_dim;      //number of dimensions
    int D_dim_0;        //number of dimensions for the original sequence
    int i, j;
    int directly_update_eta = 1; //flag to use workaround fix for errors in indirect delta_eta calculation  


    //indices used for Generate_D_Array
    int Dmax;
    int D_index;
    int codon_index;
    int aa_index;

    //double Ne = Ne;
    // double mu = MU; //mutation rate per nt
    double phi = seq -> phi_obs;        //protein production rate
    double qPhi = Q * phi;
    double two_qNePhi = 2*Ne*qPhi;
    double factor;      //ratio of mut and wt elong_pr, used to rescale sigma_vecs
    double inv_factor;  //inverse of above
    
    /* Parameters related to simulation length*/
    double evol_time=0;
    double time_step;
    double wait_parameter;
    double max_time,avg_time_step;
    int evol_steps;

    //  double sigma_ratio_vec[MAX_AA];
    double b_over_c_vec[MAX_AA];
    double delta_eta_first_term_vec[MAX_AA];
    double delta_eta_second_term_vec[MAX_AA];

    double delta_eta_vec[MAX_DDIM]; //array for storing $\Delta \eta_{i,j}$ values for all
    // one step mutants of the resident allele
    double pi_vec[MAX_DDIM];    //array for storing values of $\pi(i->j)$ for all one step mutants
    double pi_total;    //sum over pi_vec values.
    double mutation_vec[MAX_DDIM];  //array for storing values of mutation rates $\mu_{i->j}$ for all one step mutants
    double pr_vec[MAX_DDIM];    //vector of pi * mu values
    double pr_total;    //sum over pr_vec values.
    double mut_pr;      //RV representing replacement allele
    double curr_pr;     //used to find replacement allele

    //double delta_eta_list[MAX_EVOL_STEPS];
    double delta_eta_mean = 0;  //mean effect of synon sub for the current wt seq-- def differs from before
    double delta_eta_var = 0;   //var effect of synon sub for the current wt seq --definition differs from before

    double evol_delta_eta_mean = 0; //mean effect of synon sub evolution
    double evol_delta_eta_var = 0;  //var effect of synon sub evolution
    double evol_delta_mean; //used to update mean and var

    double delta_eta;   //eta_wt-eta_mut


    //variables used for updating eta and calculating delta_eta directly
    
    double previous_eta;
    double new_xi;
    double new_sigma_n;
    double new_eta;
    double direct_calc_delta_eta; //delta_eta based on direct calculations of eta for new and previous sequences.
    double abs_error_delta_eta;
    double rel_error_delta_eta;

    /*Parameters related to liklihood calculation*/
    double seq_mu_fitness_Ne = 0;
    double step_mu_fitness_Ne = 0;
    double seq_fitness_Ne = 0;
    double step_fitness_Ne = 0;
    double seq_mu = 0;
    double step_mu = 0;
    
    double step_eta;

                                                       //but i think it is right. May need to change this 2/24/13
                                                       //eta_wt - eta_mut
                                                       //So, if eta_wt > eta_mut, delta_eta_0_i is positive. 
    struct timeval curr_time;

    const gsl_rng_type *rng_type;
    gsl_rng *rn;
    

/*========================================================*
 *                1) SET UP GSL RNG                       *
 *      Create a generator chosen by the                  *
 *      environment variable GSL_RNG_TYPE                 *
 *      use default gsl for generating uniform rn         *
 *      from which other rn functions are derived         *
 *========================================================*/

    gsl_rng_env_setup();

    rng_type = gsl_rng_default;
    rn = gsl_rng_alloc(rng_type);

    //seed rng using clock
    gettimeofday(&curr_time, NULL);
    gsl_rng_set(rn, (curr_time.tv_sec + curr_time.tv_usec));    // (const gsl_rng * r, long int s)   

/*========================================================*
 *    2) CALCULATE SEQUENCE SIGMA_VEC AND B_OVER_C_VEC    *
 *========================================================*/

    Convert_Sigma_Vec_to_Sigma_Ratio_Vec(&(seq->sigma_vec[0]),
                         &(seq->sigma_ratio_vec[0]),
                         aa_count);

    Calc_B_over_C_Vec(&(seq->codon_index[0]), b_over_c_vec, aa_count);
    

    
/*========================================================*
 *              3) START SEQ ITERATIONS                   *
 * The first round of simulation steps should burn in the *
 * sequence with a fixed number of steps. The second      *
 * of simulation steps should run for a given amount of   *
 * evolution time.                                        *
 *========================================================*/

    evol_steps = 0;
    time_step = 0;
    
    //*=========*//
    //* Round 1 *//
    //*=========*//
    
    for (evol_steps = 0; evol_steps < total_evol_steps; evol_steps++) {
                
        //*=====================================================================*
        //*  4.b) Determine Dimensionality/# one step neighbors of the system.  *
        //*  D_dim DOES change over the simulation!                             *
        //*  The changes in its values are due to the effects of 6 codon aa     *
        //*  For example, for Arginine AGA can mutate to AGG or CGA.            *
        //*  In contrast, CGA can mutate to AGA, CGC, CGG, or CGU.              *
        //*=====================================================================*

        D_dim = Calc_Dimensionality(&(seq->codon_index[0]), aa_count);
                
        //*=============================================================================*
        //*  4.c) Create D_array                                                        *
        //*                                                                             *
        //*  D_array is a 2D matrix of size 2 x Dimensionality of sequence              *
        //*      Dimensionality = # one step syonymous neighbors                        *
        //*      D_array[i][0] = j, where j is the position of an amino acid in the ORF *
        //*      D_array[i][1] = k, where k is a codon_index value for an alternative   *
        //*                            synonymous codon for the wt sequence             *
        //*=============================================================================*
        
        D_index = 0;
        
        for (i = 0; i < aa_count; i++) {
            codon_index = seq->codon_index[i];
            aa_index = seq->aa_index[i];
            Dmax = Codon[codon_index].num_synonym;

            //go through each synonym
            for (j = 0; j < Dmax; j++) {

                mutation_vec[D_index] = Codon[codon_index].syn_relative_mutation_rate[j];

                (D_array[D_index][0]) = (int)(i);    //amino acid position
                (D_array[D_index++][1]) = (Codon[codon_index].synonym_index[j]);   //synonymous codon for position D_index
            }
        }     
        
		//*==========================================================*
		//*  4.d) Calculate Delta_Eta and fixation probability (pi)  *
		//*  for all one step neighbors.                             *
		//*==========================================================*
		
		Calc_Delta_Eta_NSE_First_and_Second_Term_Vecs(&(seq->sigma_ratio_vec[0]),
													  b_over_c_vec,
													  delta_eta_first_term_vec,
													  delta_eta_second_term_vec,
													  aa_count);

		Calc_Delta_Eta_NSE_Vec(&(seq->codon_index[0]), b_over_c_vec,
							   delta_eta_first_term_vec,
							   delta_eta_second_term_vec,
							   delta_eta_vec, aa_count, D_dim);

		
		Calc_Pi_Vec(delta_eta_vec, pi_vec, &pi_total, phi, qPhi,
					two_qNePhi, aa_count, D_dim);

		//*=====================================================================*
		//*  4.e) Calculate replacement probability for all one step neighbors  *
		//*       --> Test to see if there are mutation effects                 *
		//*           -Uses AT bias and transition/transversion bias OR         *
		//*            mutation rates loaded from a mutation file.              *
		//*=====================================================================*
		
		if ((mu_bias == 0.5) && gamma_ratio == 1. && 0) {
			//* flat mutation so pi values = pr values
			//* Assign pr_vec address to pi_vec address
			//* This should work so long as the dimensions of these
			//* vectors are the same
			//*       &(pr_vec[0]) = &(pi_vec[0]);//way of copying pointers according to Kernighan and Ritchie p. 98
			//* The above doesn't work because pr_vec represents a specific space in memory and it its address can't
			//* be reassigned.  It can be copied, however.

			memcpy(pr_vec, pi_vec, D_dim * sizeof(double));

			pr_total = pi_total;

		} else {
			//* Calculate true pr of replacement, weighing by mutation rates
			Calc_Replacement_Pr_Vec(mutation_vec, pi_vec, pr_vec,
						&pr_total, D_dim);
		}
       
        
        //*=============================================================*
        //*  4.f) Calculate expected time until replacement             *
        //*  We are assuming that replacement                           *
        //*  happens quickly. Ideally we would expect the wait          *
        //*  parameter to be >>1, but it likely doesn't matter          *
        //*  The time to replacement should follow a geometric dist     *
        //*  which we approximate with an exponential distribution.     *
        //*  Note the parameter for GSL's RNG argument is the expected  *
        //*  wait time, so it is 1/replacement rate                     *
        //*=============================================================*
        
        wait_parameter = 1 / (Ne * mu * pr_total);  // This is the 'failure' or leaving rate of the resident allele.
        
        //* Calculate based on wait_parameter only or random exponential
        
        if (benchmark||non_random_time_step){
            time_step = wait_parameter;
        }else{
            time_step = gsl_ran_exponential(rn, wait_parameter);
        }
       
        //*==============================================*
        //*              DEBUGGING TOOLS                 *
        //*==============================================*
        //*  Prints a detailed summary of all one step   *
        //*  neighbors.                                  *
        //*==============================================*
        
        if(1==0){
            
            Rprintf("############# Step Number: %d #############\n",evol_steps);
            Rprintf("\t\t\t\tdelta_eta\tpi_i_j\tpr_i_j\n");
            for(i=0;i<D_dim;i++){
                aa_position = D_array[i][0];
                wt_index = seq->codon_index[aa_position];
                mut_index = D_array[i][1];
                
                Rprintf("Codon %d: %f  ->  Codon %d: %f\t%g\t%g\t%g\n",wt_index,Codon[wt_index].elong_rate,mut_index,Codon[mut_index].elong_rate,delta_eta_vec[i],pi_vec[i],pr_vec[i]);
                error("Delta eta did not work\n");
            }
        }   
        
        //*==============================================*
        //*  4.g) Choose a RV to determine which allele  *
        //*  replaces the resident.                      * 
        //*==============================================*
        
        mut_pr = gsl_rng_uniform(rn)*pr_total;
        if (benchmark) mut_pr = 0.5 * pr_total;

            //* Note: pr_total is not scaled by the constants Ne or mu)
        
        //*===========================================================*
        //*  4.h) Search forward to find substitution allele.         *
        //*  This is likely more efficient than starting at the end   *
        //*  b/c there is weaker selection at the start of the        *
        //*  sequence.                                                *
        //*  Some kind of newton search could be even more efficient  *
        //*===========================================================*
        
        curr_pr = 0.0;
        i = 0;
        while (curr_pr < mut_pr) {
            curr_pr += pr_vec[i++];
        }

        i--;
        D_index = i;

		//***Store This for exchange with bridging***//
        mut_codon_pos = D_array[D_index][0];
        mut_codon_index = D_array[D_index][1];

		//*======================*
		//*  Get wt codon index  *
		//*======================*
        
        wt_codon_index = seq->codon_index[mut_codon_pos];
       
		//*==============*
		//*  Get codons  *
		//*==============*
        
        strcpy(wt_codon, Codon[wt_codon_index].codon);
        strcpy(mut_codon, Codon[mut_codon_index].codon);
        
        if(benchmark) {
            Rprintf("Substituting Codon %d %s for Codon %d %s\n",wt_codon_index,wt_codon,mut_codon_index,mut_codon);
        }
               
        //*========================================*
        //*      4.i) Update everything            *
        //*========================================*
        
		//* update codon index
        seq->codon_index[mut_codon_pos] = mut_codon_index;

		//* update codon counts
        seq->codon_cts[wt_codon_index]--;
        seq->codon_cts[mut_codon_index]++;

		//* FIXME: This delta_eta appears to be slightly off!
		//* Using temporary work around
        delta_eta = delta_eta_vec[i];   //should be \eta_i - \eta_j

        if(directly_update_eta){
            //* WORKAROUND: Calculate eta for new sequence from scratch
            //* this line overwrites the updating of sigma_vec done above which might be the cause of our errors in delta_eta
            //* it should be possible to tell since the delta_eta would work for the first evolutionary step but not the later ones.
            //* TODO: Comment out this line to see if there's a problem with how we update sigma_vec above
            
            new_eta = Calc_Eta_NSE (seq);
            //* this overwrites update of seq->eta_obs using delta_eta
            seq->eta_obs = new_eta;

            if (1 == 1) {
              //* check to make sure delta_eta is working right
              direct_calc_delta_eta = (previous_eta - new_eta);
              abs_error_delta_eta = abs(direct_calc_delta_eta - delta_eta);
              rel_error_delta_eta = abs_error_delta_eta/direct_calc_delta_eta;
              
              //Rprintf("Step %d: Position: %d Codon %d %s --> Codon %d %s\td_eta: %f \tdir_d_eta: %f\n",evol_steps,mut_codon_pos,wt_codon_index,wt_codon,mut_codon_index,mut_codon,delta_eta,direct_calc_delta_eta);
              if (rel_error_delta_eta > 1.0E-6) {           
                 Rprintf ("ERROR: direct_d_eta: %f != delta_eta: %f \n\tAA: %c Codon %d %s --> Codon %d %s\n", direct_calc_delta_eta, delta_eta,Codon[wt_codon_index].aa,wt_codon_index,Codon[wt_codon_index].codon,mut_codon_index,Codon[mut_codon_index].codon);
                 Rprintf("############# Step Number: %d D_dim: %d #############\n",evol_steps,D_dim);
			     Rprintf("\t\t\t\tdelta_eta\tpi_i_j\tpr_i_j\n");
				 for(i=0;i<D_dim;i++){
					aa_position = D_array[i][0];
					wt_index = seq->codon_index[aa_position];
					mut_index = D_array[i][1];
                
					Rprintf("Codon %d: %f  ->  Codon %d: %f\t%g\t%g\t%g\n",wt_index,Codon[wt_index].elong_rate,mut_index,Codon[mut_index].elong_rate,delta_eta_vec[i],pi_vec[i],pr_vec[i]);
					
				 }
				 error("Delta eta did not work\n");
              }
            }
            //Rprintf("delta_eta: %f\tdirect_delta_eta: %f\n",delta_eta,(previous_eta - new_eta));
            delta_eta = (previous_eta - new_eta);

        }else{

              //* Update Sigma_Ratio_Vec by dividing the values up to  
              //* mut_codon_pos by elong_pr_i/elong_pr_j               
              //* values. After this point don't change since they     
              //* involve the factor in both the                       
              //* denominator and numerator.                           
          
          factor = Codon[mut_codon_index].elong_pr /
            Codon[wt_codon_index].elong_pr;
          inv_factor =
            Codon[wt_codon_index].elong_pr /
            Codon[mut_codon_index].elong_pr;

          for (i = 1; i <= mut_codon_pos; i++) {
            seq->sigma_ratio_vec[i] *= inv_factor;
          }
          
              //* update sigma_vec for values at and above mut_codon_pos;
          for (i = mut_codon_pos; i < aa_count; i++) {          
            seq->sigma_vec[i] *= factor;
          }
          
              //* update eta_obs by subtracting difference: \eta_wt - (\eta_wt - \eta_mut) = \eta_mut
          seq->eta_obs -= delta_eta;

        } //end else


		//* update b_over_c_vec
        b_over_c_vec[mut_codon_pos] = B / Codon[mut_codon_index].elong_rate;

        //* update evol_delta_mean/var vals
        //evol_delta_mean = delta_eta - delta_eta_mean;
        //evol_delta_eta_mean += evol_delta_mean / (evol_time + time_step);
        //evol_delta_eta_var +=  evol_delta_mean * (delta_eta - evol_delta_eta_mean);

		//* update values for printing
        //seq->evol_delta_eta_mean = evol_delta_eta_mean;
        //seq->evol_delta_eta_var = evol_delta_eta_var;
        
        evol_time += time_step; //increment evol_time (changed for loop to maximum#steps)
    }
    
    //* ROUND 2
    if (global_max_time < 0) {
        avg_time_step = evol_time/total_evol_steps;//average time spent at each step
		max_time = -avg_time_step*aa_count*global_max_time;
		
    } else {
        max_time = global_max_time * avg_time_step;
    }
    
    
    for (evol_time = 0; evol_time < max_time; evol_time += time_step) { 
        evol_steps++; 
               
        //*=====================================================================*
        //*  4.b) Determine Dimensionality/# one step neighbors of the system.  *
        //*  D_dim DOES change over the simulation!                             *
        //*  The changes in its values are due to the effects of 6 codon aa     *
        //*  For example, for Arginine AGA can mutate to AGG or CGA.            *
        //*  In contrast, CGA can mutate to AGA, CGC, CGG, or CGU.              *
        //*=====================================================================*

        D_dim = Calc_Dimensionality(&(seq->codon_index[0]), aa_count);
        
        if (evol_steps == 0) {
            D_dim_0 = D_dim;    //store original dimensionality as size of synonymous space for estimation of denominator
            //Rprintf("Seq ID: %s\tD_dim_0:%d\n",seq->id,D_dim_0);
        }
        
        //*=============================================================================*
        //*  4.c) Create D_array                                                        *
        //*                                                                             *
        //*  D_array is a 2D matrix of size 2 x Dimensionality of sequence              *
        //*      Dimensionality = # one step syonymous neighbors                        *
        //*      D_array[i][0] = j, where j is the position of an amino acid in the ORF *
        //*      D_array[i][1] = k, where k is a codon_index value for an alternative   *
        //*                            synonymous codon for the wt sequence             *
        //*=============================================================================*
        
        D_index = 0;
        for (i = 0; i < aa_count; i++) {
            codon_index = seq->codon_index[i];
            aa_index = seq->aa_index[i];
            Dmax = Codon[codon_index].num_synonym;

            //go through each synonym

            for (j = 0; j < Dmax; j++) {

                mutation_vec[D_index] =
                    Codon[codon_index].
                    syn_relative_mutation_rate[j];

                (D_array[D_index][0]) = (int)(i);    //amino acid position
                (D_array[D_index++][1]) = (Codon[codon_index].synonym_index[j]);   //synonymous codon for position D_index
            }
        }
        
        //*==========================================================*
		//*  4.d) Calculate Delta_Eta and fixation probability (pi)  *
		//*  for all one step neighbors.                             *
		//*==========================================================*
		
		Calc_Delta_Eta_NSE_First_and_Second_Term_Vecs(&(seq->sigma_ratio_vec[0]),
													  b_over_c_vec,
													  delta_eta_first_term_vec,
													  delta_eta_second_term_vec,
													  aa_count);

		Calc_Delta_Eta_NSE_Vec(&(seq->codon_index[0]), b_over_c_vec,
							   delta_eta_first_term_vec,
							   delta_eta_second_term_vec,
							   delta_eta_vec, aa_count, D_dim);

		
		Calc_Pi_Vec(delta_eta_vec, pi_vec, &pi_total, phi, qPhi,
					two_qNePhi, aa_count, D_dim);

		//*=====================================================================*
		//*  4.e) Calculate replacement probability for all one step neighbors  *
		//*       --> Test to see if there are mutation effects                 *
		//*           -Uses AT bias and transition/transversion bias OR         *
		//*            mutation rates loaded from a mutation file.              *
		//*=====================================================================*
		
		if ((mu_bias == 0.5) && gamma_ratio == 1. && 0) {
			//* flat mutation so pi values = pr values
			//* Assign pr_vec address to pi_vec address
			//* This should work so long as the dimensions of these
			//* vectors are the same
			//*       &(pr_vec[0]) = &(pi_vec[0]);//way of copying pointers according to Kernighan and Ritchie p. 98
			//* The above doesn't work because pr_vec represents a specific space in memory and it its address can't
			//* be reassigned.  It can be copied, however.

			memcpy(pr_vec, pi_vec, D_dim * sizeof(double));

			pr_total = pi_total;

		} else {
			//* Calculate true pr of replacement, weighing by mutation rates
			Calc_Replacement_Pr_Vec(mutation_vec, pi_vec, pr_vec,
						&pr_total, D_dim);
		}
       
        
        //*=============================================================*
        //*  4.f) Calculate expected time until replacement             *
        //*  We are assuming that replacement                           *
        //*  happens quickly. Ideally we would expect the wait          *
        //*  parameter to be >>1, but it likely doesn't matter          *
        //*  The time to replacement should follow a geometric dist     *
        //*  which we approximate with an exponential distribution.     *
        //*  Note the parameter for GSL's RNG argument is the expected  *
        //*  wait time, so it is 1/replacement rate                     *
        //*=============================================================*
        
        wait_parameter = 1 / (Ne * mu * pr_total);  // This is the 'failure' or leaving rate of the resident allele.
        
        //* Calculate based on wait_parameter only or random exponential
        
        if (benchmark||non_random_time_step){
            time_step = wait_parameter;
        }else{
                time_step = gsl_ran_exponential(rn, wait_parameter);
        }
        
        //*==============================================*
        //*              DEBUGGING TOOLS                 *
        //*==============================================*
        //*  Prints a detailed summary of all one step   *
        //*  neighbors.                                  *
        //*==============================================*
        
        if(1==0){
            
            Rprintf("############# Step Number: %d #############\n",evol_steps);
            Rprintf("\t\t\t\tdelta_eta\tpi_i_j\tpr_i_j\n");
            for(i=0;i<D_dim;i++){
                aa_position = D_array[i][0];
                wt_index = seq->codon_index[aa_position];
                mut_index = D_array[i][1];
                
                Rprintf("Codon %d: %f  ->  Codon %d: %f\t%g\t%g\t%g\n",wt_index,Codon[wt_index].elong_rate,mut_index,Codon[mut_index].elong_rate,delta_eta_vec[i],pi_vec[i],pr_vec[i]);
                error("Delta eta did not work\n");
            }
        }   
        
        //*==============================================*
        //*  4.g) Choose a RV to determine which allele  *
        //*  replaces the resident.                      * 
        //*==============================================*
        
        mut_pr = gsl_rng_uniform(rn)*pr_total;
        if (benchmark) mut_pr = 0.5 * pr_total;

            //* Note: pr_total is not scaled by the constants Ne or mu)
        
        //*===========================================================*
        //*  4.h) Search forward to find substitution allele.         *
        //*  This is likely more efficient than starting at the end   *
        //*  b/c there is weaker selection at the start of the        *
        //*  sequence.                                                *
        //*  Some kind of newton search could be even more efficient  *
        //*===========================================================*
        
        curr_pr = 0.0;
        i = 0;
        while (curr_pr < mut_pr) {
            curr_pr += pr_vec[i++];
        }

        i--;
        D_index = i;

		//***Store This for exchange with bridging***//
        mut_codon_pos = D_array[D_index][0];
        mut_codon_index = D_array[D_index][1];

		//*======================*
		//*  Get wt codon index  *
		//*======================*
        
        wt_codon_index = seq->codon_index[mut_codon_pos];
       
		//*==============*
		//*  Get codons  *
		//*==============*
        
        strcpy(wt_codon, Codon[wt_codon_index].codon);
        strcpy(mut_codon, Codon[mut_codon_index].codon);
        
        if(benchmark) {
            Rprintf("Substituting Codon %d %s for Codon %d %s\n",wt_codon_index,wt_codon,mut_codon_index,mut_codon);
        }
               
        //*========================================*
        //*      4.i) Update everything            *
        //*========================================*
        
		//* update codon index
        seq->codon_index[mut_codon_pos] = mut_codon_index;

		//* update codon counts
        seq->codon_cts[wt_codon_index]--;
        seq->codon_cts[mut_codon_index]++;

		//* FIXME: This delta_eta appears to be slightly off!
		//* Using temporary work around
        delta_eta = delta_eta_vec[i];   //should be \eta_i - \eta_j

        if(directly_update_eta){
            //* WORKAROUND: Calculate eta for new sequence from scratch
            //* this line overwrites the updating of sigma_vec done above which might be the cause of our errors in delta_eta
            //* it should be possible to tell since the delta_eta would work for the first evolutionary step but not the later ones.
            //* TODO: Comment out this line to see if there's a problem with how we update sigma_vec above
            
            new_eta = Calc_Eta_NSE (seq);
            //* this overwrites update of seq->eta_obs using delta_eta
            seq->eta_obs = new_eta;

            if (1 == 1) {
              //* check to make sure delta_eta is working right
              direct_calc_delta_eta = (previous_eta - new_eta);
              abs_error_delta_eta = abs(direct_calc_delta_eta - delta_eta);
              rel_error_delta_eta = abs_error_delta_eta/direct_calc_delta_eta;
              
              //Rprintf("Step %d: Position: %d Codon %d %s --> Codon %d %s\td_eta: %f \tdir_d_eta: %f\n",evol_steps,mut_codon_pos,wt_codon_index,wt_codon,mut_codon_index,mut_codon,delta_eta,direct_calc_delta_eta);
              if (rel_error_delta_eta > 1.0E-6) {           
                 Rprintf ("ERROR: direct_d_eta: %f != delta_eta: %f \n\tAA: %c Codon %d %s --> Codon %d %s\n", direct_calc_delta_eta, delta_eta,Codon[wt_codon_index].aa,wt_codon_index,Codon[wt_codon_index].codon,mut_codon_index,Codon[mut_codon_index].codon);
                 Rprintf("############# Step Number: %d D_dim: %d #############\n",evol_steps,D_dim);
			     Rprintf("\t\t\t\tdelta_eta\tpi_i_j\tpr_i_j\n");
				 for(i=0;i<D_dim;i++){
					aa_position = D_array[i][0];
					wt_index = seq->codon_index[aa_position];
					mut_index = D_array[i][1];
                
					Rprintf("Codon %d: %f  ->  Codon %d: %f\t%g\t%g\t%g\n",wt_index,Codon[wt_index].elong_rate,mut_index,Codon[mut_index].elong_rate,delta_eta_vec[i],pi_vec[i],pr_vec[i]);
					
				 }
				 error("Delta eta did not work\n");
              }
            }
            //Rprintf("delta_eta: %f\tdirect_delta_eta: %f\n",delta_eta,(previous_eta - new_eta));
            delta_eta = (previous_eta - new_eta);

        }else{

		  //* Update Sigma_Ratio_Vec by dividing the values up to  
		  //* mut_codon_pos by elong_pr_i/elong_pr_j               
		  //* values. After this point don't change since they     
		  //* involve the factor in both the                       
		  //* denominator and numerator.                           
	  
          factor = Codon[mut_codon_index].elong_pr /
            Codon[wt_codon_index].elong_pr;
          inv_factor =
            Codon[wt_codon_index].elong_pr /
            Codon[mut_codon_index].elong_pr;

          for (i = 1; i <= mut_codon_pos; i++) {
            seq->sigma_ratio_vec[i] *= inv_factor;
          }
          
              //* update sigma_vec for values at and above mut_codon_pos;
          for (i = mut_codon_pos; i < aa_count; i++) {          
            seq->sigma_vec[i] *= factor;
          }
          
              //* update eta_obs by subtracting difference: \eta_wt - (\eta_wt - \eta_mut) = \eta_mut
          seq->eta_obs -= delta_eta;

        } //end else


		//* update b_over_c_vec
        b_over_c_vec[mut_codon_pos] = B / Codon[mut_codon_index].elong_rate;

        //* update evol_delta_mean/var vals
        //evol_delta_mean = delta_eta - delta_eta_mean;
        //evol_delta_eta_mean += evol_delta_mean / (evol_time + time_step);
        //evol_delta_eta_var +=  evol_delta_mean * (delta_eta - evol_delta_eta_mean);

		//* update values for printing
        //seq->evol_delta_eta_mean = evol_delta_eta_mean;
        //seq->evol_delta_eta_var = evol_delta_eta_var;

        
    }

/*========================================================*
 *                 END EVOLUTION ITERATIONS               *
 *========================================================*/
    
    

    //update Seq[]
    seq->evol_steps = evol_steps;
    seq->evol_time = evol_time;
    seq->delta_eta_mean = delta_eta_mean;
    seq->delta_eta_var = delta_eta_var;

    //Finish fixation pr demoninator approximation
    seq_mu /= evol_time;
    seq_mu_fitness_Ne /= evol_time;
    seq_fitness_Ne /= evol_time;

    gsl_rng_free(rn);

    return 1;
}

void MCMC_Sequence(struct seq_struct *seq){
        
	int i;
	int recent_steps[250]; 	//this is a running buffer of acceptances and rejections. 
							//it is indexed each step as recent_steps[step_num%250] so
							//that each successive step overwrites the step 250 steps ago.
	
	int accept_num=0,step_num=0;
	double accept_ratio;
	double desired_acc_ratio = 0.3;
	double accept_pr;
	
	struct seq_struct tSeq;
	
			
	int num_sub=5;		
	int max_sub = 50;
	double randNums[200];
	double step_eta,step_mu;
	
	struct timeval curr_time;
	
	const gsl_rng_type *rng_type;
    gsl_rng *rn;
    		
	/*========================================================*
	 *                1) SET UP GSL RNG                       *
	 *      Create a generator chosen by the                  *
	 *	    environment variable GSL_RNG_TYPE                 *
	 *      use default gsl for generating uniform rn         *
	 *      from which other rn functions are derived         *
	 *========================================================*/

    gsl_rng_env_setup();

    rng_type = gsl_rng_default;
    rn = gsl_rng_alloc(rng_type);

    //seed rng using clock
    gettimeofday(&curr_time, NULL);
    gsl_rng_set(rn, (curr_time.tv_sec + curr_time.tv_usec));    // (const gsl_rng * r, long int s)   		
        
    
    ///////////////////////////////////
    //1: Initialize tSeq.codon_index //
    ///////////////////////////////////
    
    for(i=0;i<seq->aa_count;i++){
		tSeq.codon_index[i] = seq->codon_index[i];
	}
	tSeq.aa_count = seq->aa_count;
    
    
    
    /////////////////////////////////////////////////////
    //2: Start for loop-- start num substitutions at 5 //
    /////////////////////////////////////////////////////
    
	while(step_num < total_evol_steps) {
		
	//////////////////////////////////	
    //3: Take a step in codon space //
    //////////////////////////////////
    
		for(i=0;i<2*num_sub;i++){
			*(randNums + i) = gsl_rng_uniform(rn);
		}
		for(i=0;i<seq->aa_count;i++){
			tSeq.codon_index[i] = seq->codon_index[i];
		}
		
		flip_codons(seq,&(tSeq),randNums,num_sub);
    
    ///////////////////////////
    //4: Calc new eta and mu //
    ///////////////////////////
		
		tSeq.eta_initial = tSeq.eta_obs = Calc_Eta_NSE(&tSeq);
		
		Codon_Counts(tSeq.codon_index,tSeq.codon_cts,tSeq.aa_count);
		
		tSeq.mu_obs = Calc_Seq_Mu(tSeq.codon_cts,tSeq.aa_count);
	
	////////////////////////////////////////	
    //5: accept pr = min(w(new)/w(old),1) //
    ////////////////////////////////////////
    
		accept_pr = fmin(tSeq.mu_obs/seq->mu_obs*exp(-Q*Ne*seq->phi_obs*(tSeq.eta_obs-seq->eta_obs)),1);
		//Rprintf("num_sub: %d curr_eta: %f new_eta: %f accept_pr: %f \n",num_sub,seq->eta_obs,tSeq.eta_obs,accept_pr);
    
    /////////////////////////////////
    //7: Draw rn and accept/reject //
    /////////////////////////////////
    
		*randNums = gsl_rng_uniform(rn);
		
		if(*randNums < accept_pr){
			
			step_eta = tSeq.eta_obs;
			step_mu = tSeq.mu_obs;
			
			
			//Replace seq with tSeq - only need to change codon_cts, codon_index, eta_obs, mu_obs
			seq->eta_obs = step_eta;
			seq-> mu_obs = step_mu;
			for(i=0;i<64;i++){
				*(seq->codon_cts+i) = tSeq.codon_cts[i];
			}
			for(i=0;i<seq->aa_count;i++){
				*(seq->codon_index+i) = tSeq.codon_index[i];
			}
			recent_steps[step_num%250] = 1; //keeps a running buffer of accepts and rejects
			step_num++;
			
		}else{ //If this is not a new step, only increment time
			recent_steps[step_num%250] = 0;
			step_num++;
		}
    
    ////////////////////////////////////////////////////////////////////////////
    //8: adjust num_substitutions											  //
    //      -if(accept ratio > desired_acc_ratio): increase num substitutions //
    //      -if(accept ratio < desired_acc_ratio): decrease num substitutions //
    ////////////////////////////////////////////////////////////////////////////
    
		if(step_num > 250 && step_num%50==0){
			accept_num = 0;
			//tally up all the acceptances from last 250 steps
			for(i=0;i<250;i++){
				accept_num+=recent_steps[i];
			}
			
			accept_ratio = (double) accept_num/250;
			
			if((accept_ratio > (desired_acc_ratio+0.05))&&(num_sub<max_sub)){
				num_sub++;
			}else if((accept_ratio < (desired_acc_ratio-0.05))&&(num_sub>1)){
				num_sub--;
			}
		}
		//Rprintf("step_eta: %f ",seq->eta_obs);
		//Rprintf("accept_num: %d step_num: %d ",accept_num,step_num);
		//Rprintf("accept_ratio: %f ",accept_ratio);
		//Rprintf("num_sub: %d\n",num_sub);
		
    
    
	}
    
    //Rprintf("num_sub: %d\n",num_sub);
    
	gsl_rng_free(rn);
	
}



void CES_new(int *CODON_INDEX,                  double *ELONGATION_RATES,         double *MUTATION_RATES, 
			 int *AA_COUNT,                     double *PHI,                      double *POP,
			 double *A_1,                       double *A_2,                      double *AT_BIAS,                  
			 int *BENCH,                        double *BEE,                      double *GMT,                      
			 int *IGNORE,                       double *GAMMA,                    double *SCALE_FACTOR, 
			 int *MES,             
			 int *BIS,                          char ** AA_VEC,                   char ** CODON_VEC,
			 int * N_CODON,						char ** SIMULATION_METHOD)	
{
	
    int i, j;
    int aa_count;
    double max_time;
    double eta_original;

    
    //* Initializing the random number generator
    srand(time(NULL));  //should this be NULL? mikeg
    

    //*============================*
    //*  1) Process input from R   *
    //*============================*

    //* Import Sequence data and genome parameters
    Process_R_input(CODON_INDEX,ELONGATION_RATES,MUTATION_RATES,AA_COUNT,
                    PHI, POP, A_1, A_2, AT_BIAS,
                    BEE, IGNORE, GAMMA, SCALE_FACTOR,
                    AA_VEC, CODON_VEC, N_CODON);
    
    //* Simulation parameters
    benchmark = *BENCH;
    global_max_time = *GMT;
    simulation_steps = *MES;
    burn_in_steps = (*BIS) * Sequence.aa_count;
    
    
    Process_AA_Information();

    //* initialize codon structures
    Generate_Codon_Structures();
    
    //* adjust aa_counts if we are ignoring aa's on tail end
    //* the idea is that a protein can still function if it is missing
    //* the last few amino acids (i.e. translation has a nonsense error
    //* on the last few codons). We do not use this in the default setting
    if (ignore_aa > 0) {
        for (i = 0; i < n_seq; i++) {
            Sequence.aa_count -= ignore_aa;
        }
    }
    aa_count = Sequence.aa_count;
    
   
    //*=========================================================*
    //*  2) Calculate length of simulation                      *
    //*     --Default behavior is to use fixed number of        *
    //*       steps, calculated from *MES and *BIS, but we can  *
    //*       also used fixed amount of time, calculated from   *
    //*       *GMT                                              *
    //*=========================================================*

    
    total_evol_steps = burn_in_steps + simulation_steps; //BIS and MES parameters passed from R
        
    Sequence.max_time = max_time;

    //*==================================*
    //*  3) Calculate initial eta and mu *
    //*==================================*
    
    Sequence.eta_initial = Sequence.eta_obs =
        Calc_Eta_NSE(&Sequence);
        
    Codon_Counts(Sequence.codon_index,Sequence.codon_cts,Sequence.aa_count);
		
	Sequence.mu_obs = Calc_Seq_Mu(Sequence.codon_cts,Sequence.aa_count);

    //*========================*
    //* 4) Simulate Evolution  *
    //*========================*
    
    switch(**SIMULATION_METHOD){
    case 'E':
		Evolve_Sequence(&Sequence);
		break;
		
	case 'M':
		MCMC_Sequence(&Sequence);
		break;
	}
	
    //*=========================*
    //* 5) Update *CODON_INDEX  *
    //*=========================*
    
    for(j=0;j<*AA_COUNT;j++) {
        *(CODON_INDEX + j) = Sequence.codon_index[j];
    }
            
}


double Calc_Eta_NSE_R_wrapper (int * CODON_INDEX, int * AA_COUNT, double * ELONG_PR, int * N_CODONS, double * A_1, double * A_2,double * ETA){
	double xi, eta;
	int i;
	
	//1) Set up sequence structure
	//1.a) codon_index
	for(i=0;i<*AA_COUNT;i++){
		Sequence.codon_index[i] = CODON_INDEX[i];
	}
	//1.b) aa_count
	Sequence.aa_count = *AA_COUNT;
	
	//2) Set up codon elong_pr's
	for(i=0;i<*N_CODONS;i++){
		Codon[i].elong_pr = ELONG_PR[i];
	}
	
	//3) Set up initiation and elongation cost 
	A1 = *A_1;
	A2 = *A_2;
	
	*ETA = Calc_Eta_NSE(&Sequence);
	
	
	
}

/*========================================*
 *      CALCULATE ENTIRE DENOMINATOR      *
 * This function will allow us to compare *
 * our estimates of the lik_denom to real *
 * value.                                 *
 *========================================*/
 
 void calc_exact(int *CODON_INDEX,                  int *AA_COUNT,                    double *PHI,
                 double *ELONGATION_RATES,          double *MUTATION_RATES,           double *POP,                      
                 double *A_1,                       double *A_2,                      double *AT_BIAS,
                 double *BEE,                       int *IGNORE,                      double *GAMMA,                    
                 double *SCALE_FACTOR,              double *LIK,                      double *ETA_MEAN,
                 double *ETA_VAR,                   double *ETA_OBS,                  double *ETA_MIN,
                 double *ETA_MAX,                   double *BINS,                     int *NUM_BINS,
                 char **AA_VEC,                      char **CODON_VEC,                 int *N_CODON){

    
    int i;
    int increment_next=0, counter=0;
    int aa_index,aa_position,aa_count;
    int num_codons;
    int size_codon_space = 1;
    int t_codon_cts[61] = {0} ,t_codon_cts_rel[61] = {0};
    int t_codon_index[MAX_AA]; //holds temporary codon_index
    
    double t_sigma_n,t_xi,t_eta,t_mu,t_fitness,t_mu_fitness; 
    double fitness_total = 0, mu_fitness_total=0, eta_total=0, eta_sq_total=0;
    double t_sigma_vec[MAX_AA];
    
    double pcnt=0;
    
    double answer[2];
    
    double * bin_lims;
    bin_lims = malloc(sizeof(double)*(*(NUM_BINS)+1));

    //*=================================*
    //* 1) Process input passed from R  *
    //*=================================*
    
    Process_R_input(CODON_INDEX,ELONGATION_RATES,MUTATION_RATES,AA_COUNT,
                        PHI, POP, A_1, A_2, AT_BIAS,
                        BEE, IGNORE, GAMMA, SCALE_FACTOR, AA_VEC,
                        CODON_VEC, N_CODON);

    Process_AA_Information();

    //* initialize codon structures
    Generate_Codon_Structures();
    
    aa_count=Sequence.aa_count;
    
    //*=================================*
    //* 2) Find size of codon space     *
    //*=================================*

    Convert_Codon_Index_to_AA_Index(Sequence.codon_index,
                            Sequence.aa_index, aa_count);

    for(i=0;i<aa_count;i++){
        aa_index = Sequence.aa_index[i];
        num_codons = AA[aa_index].num_codons;
        size_codon_space *= num_codons;
    }
    
    //*=====================================*
    //* 3) Find eta min/max and define bin  *
    //*    limits for eta hist              *
    //*=====================================*
        
    eta_min_max(Sequence.aa_index,Sequence.aa_count,answer);
    *ETA_MIN = *(answer);
    *ETA_MAX = *(answer + 1);

    intervals(*ETA_MIN,*ETA_MAX,*(NUM_BINS)+1,bin_lims);
    
    //*=====================================*
    //* 4) Initialize observed sequence and *
    //*    create temporary codon index     *
    //*=====================================*
        
        //* 4.1) Calculate eta_obs
    Sequence.eta_obs = Calc_Eta_NSE(&Sequence);
    
    //Rprintf("Eta_obs: %f\tSigma_obs: %f\txi_obs: %f\taa_count: %d\n",Sequence.eta_obs,Sequence.sigma_obs,Sequence.xi_obs,aa_count);
        
        //* 4.2) Find initial codon counts
    Codon_Counts(Sequence.codon_index,Sequence.codon_cts,aa_count);
    
    //Rprintf("Sequence.sigma_obs: %f\nSequence.xi_obs %f\nSequence.eta_obs %f\n",Sequence.sigma_obs,Sequence.xi_obs,Sequence.eta_obs);
    
        //* 4.3) Initialize temporary codon index (this will be incremented in the for loop to 
        //*      represent each permutation of the codon sequence). 
    for(i=0;i<aa_count;i++){
        aa_index = Sequence.aa_index[i];
        t_codon_index[i] = AA[aa_index].codon_index[0];
    }
    
        //* 4.4) Find codon counts of t_codon_index
    Codon_Counts(t_codon_index,t_codon_cts,aa_count);
    
        //* 4.5) Find (t_codon_cts - Sequence.codon_cts): This will give the codon counts relative
        //*      to the original seqeunce. When we calculate Mu from these counts, it will give us
        //*      mu_sim/mu_obs. 
    for(i=0;i<61;i++){ //assume no stop codons
        t_codon_cts_rel[i] = t_codon_cts[i] - Sequence.codon_cts[i];
    }
    
    //*=====================================*
    //* 5) Big ass for loop to calculate    *
    //*    fitness of every codon sequence  *
    //*    and keep running total           *
    //*=====================================*
    
    //Rprintf("Calculating denominator with %d synonymous sequences...\n",size_codon_space);
    

    for(i=0;i<size_codon_space;i++) {
        
        //*5.1) Calc Eta

          t_eta = Calc_Eta_NSE (&Sequence);

          
        //* 5.3) Calculate Mu*Fitness
        //*      Fitness = exp(-Q*Ne*phi*eta) or exp(-y*eta)
        //*      Mu=\prod_{i=1}^n(\mu_i)
        
          t_mu = Calc_Seq_Mu(t_codon_cts_rel,aa_count); //This gives mu_sim/mu_obs
        
          t_fitness = exp(-Q*Ne*Sequence.phi_obs*(t_eta-Sequence.eta_obs)); //calculate relative to eta_obs
                                                                            //to avoid very small number                                                                 
          t_mu_fitness = t_mu*t_fitness;
          
          //* Bin eta
          bin_eta(t_eta,BINS,bin_lims,t_mu_fitness,*NUM_BINS);
          
        //* 5.4) Add fitness to total fitness
        //*      Denominator of liklihood function is the sum of all 
        //*      fitnesses in codon space
        //*      Also add eta to total eta to find avg eta in codon space
          
          fitness_total += t_fitness;
          mu_fitness_total += t_mu_fitness;
          eta_total += t_eta;
          eta_sq_total += t_eta*t_eta;
        
        
        //* 5.5) Increment codon sequence
        //*      This algorithm increments the codon sequence by incrementing the codon index
        //*      of the first position. If the first position is at the highest codon index
        //*      for that amino acid, it resets the first position to the lowest codon index
        //*      for that amino acid and increments the second position. If the second position 
        //*      is at the highest codon index for that amino acid... and so on. This works 
        //*      similar to an old-fashioned rolling wheel counter.
        //*      We also need to update codon counts here...
        
          aa_position=0;
          
          do {
              increment_next=0; //assume that we do not increment the next position
              aa_index = Sequence.aa_index[aa_position];
              num_codons = AA[aa_index].num_codons;
              if(t_codon_index[aa_position] < AA[aa_index].codon_index[num_codons-1]){
              t_codon_cts_rel[t_codon_index[aa_position]++]--; //This line subtracts one from the codon count of the current codon
                                                               //and adds one to the current codon index
              t_codon_cts_rel[t_codon_index[aa_position]]++;   //This line adds one to the codon count of the next codon
              }else{
                    t_codon_cts_rel[t_codon_index[aa_position]]--; //Subtract 1 from codon count of current codon
                    t_codon_index[aa_position] = AA[aa_index].codon_index[0];
                    t_codon_cts_rel[t_codon_index[aa_position]]++; //Add 1 to codon count of next codon
                    increment_next = 1; //turn on increment flag 
              }
              aa_position++; //
          }while(increment_next && aa_position < aa_count);

          counter++;
          
          
          if(counter>size_codon_space/(pcnt*(double)size_codon_space/100)){
              //Rprintf('%d Percent Complete.\n',pcnt);
              pcnt++;
          }
    }
    
    
    //*==========================*
    //* 6) Calculate likelihood  *
    //*==========================*
         
    double likelihood,eta_mean,eta_var;
     
    likelihood = 1/(mu_fitness_total/size_codon_space);
    eta_mean = eta_total/size_codon_space;
    eta_var =  eta_sq_total/size_codon_space - pow(eta_mean,2); //E[x^2] - E[x]^2
          
    *LIK = likelihood;
    *ETA_MEAN = eta_mean;
    *ETA_VAR = eta_var;
    *ETA_OBS = Sequence.eta_obs;
    
    free(bin_lims);
}

void calc_moments(int *CODON_INDEX,                  int *AA_COUNT,                    double *PHI,
                  double *ELONGATION_RATES,          double *MUTATION_RATES,           double *POP,                      
                  double *A_1,                       double *A_2,                      double *AT_BIAS,
                  double *BEE,                       int *IGNORE,                      double *GAMMA,                    
                  double *SCALE_FACTOR,              double *ETA_MEAN,                 double *ETA_VAR, 
                  double *ETA_MIN,                   double *ETA_MAX,                  double *ETA_OBS,
                  char ** AA_VEC,                     char ** CODON_VEC,                int * N_CODON){
    //*
    //* Purpose: Calculate the moments (mean and variance) of the eta distribution analytically
    //* Math: (Put this directly into latex to view)
    //*
    //*     \documentclass{article}
    //*     \usepackage{fullpage}
    //*     \begin{document}
    
    //*     $\displaystyle{E[\eta(\vec{c})]=b\sum_{i=1}^n(a_1+a_2i)E\left[\frac{1}{c_i}\right]\prod_{j=i+1}^nE\left[\frac{c_j+b}{c_j}\right]}$

    //*     $\displaystyle{Var(\eta(\vec{c})) = b^2\left[\sum_{i=1}^n(a_1+a_2i)^2Var(Y_i) + 2\sum_{i=1}^n(a_1+a_2i)\sum_{j=i+1}^n(a_1+a_2j)Cov(Y_i,Y_j)\right]}$

    //*     $\displaystyle{Var(Y_i) = E\left[\left(\frac{1}{c_i}\right)^2\right]\prod_{j=i+1}^nE\left[\left(\frac{c_j+b}{c_j}\right)^2\right]-\left(E\left[\frac{1}{c_i}\right]\prod_{j=i+1}^nE\left[\frac{c_j+b}{c_j}\right]\right)^2}$

    //*     $\displaystyle{Cov(Y_i,Y_j)=E[Y_i*Y_j]-E[Y_i]E[Y_j]}$

    //*     $\displaystyle{Cov(Y_i,Y_j)=E[Y_i]\left(\frac{1}{\left(\prod_{k=j}^nE\left[\prod{c_k+b}{c_k}\right]\right)\left(E\left[\frac{1}{c_j}\frac{c_j+b}{c_j}\right]\right)\left(\prod_{k=j+1}^nE\left[\left(\prod_{k=j+1}^n{c_k+b}{c_k}\right)^2\right]\right)}-E(Y_j)\right)}$

    //*     $\displaystyle{H_j=\left(\prod_{k=j}^nE\left[\prod{c_k+b}{c_k}\right]\right)\left(E\left[\frac{1}{c_j}\frac{c_j+b}{c_j}\right]\right)\left(\prod_{k=j+1}^nE\left[\left(\prod_{k=j+1}^n{c_k+b}{c_k}\right)^2\right]\right)}$

    //*     \end{document}  
    
    int i;
    int aa_count,aa_index_i,aa_index_ip1;
    double prod_e[MAX_AA]={0}, sq_prod_e[MAX_AA]={0},prod_e_sq[MAX_AA]={0}; // prod_e[i] := \prod_{k=i+1}^n\frac{c_k+b}{c_k}
    double Y[MAX_AA]={0},varY[MAX_AA]={0},H[MAX_AA]={0};
    double inner_loop_sum;
    double var_eta,mean_eta=0;
    double a1pa2i;
    double answer[2];


    //*=================================*
    //* 1) Process input passed from R  *
    //*=================================*
    
    Process_R_input(CODON_INDEX,ELONGATION_RATES,MUTATION_RATES,AA_COUNT,
                        PHI, POP, A_1, A_2, AT_BIAS,
                        BEE, IGNORE, GAMMA, SCALE_FACTOR,AA_VEC,
                        CODON_VEC,N_CODON);


    Process_AA_Information(); //NOTE: this function call now also includes Generate_Codon_Structures()
    
    aa_count=Sequence.aa_count;
        
   
    //* Only need AA sequence now since overall eta distribution is independent of 
    //* observed sequence
    Convert_Codon_Index_to_AA_Index(Sequence.codon_index,
                            Sequence.aa_index, aa_count);
    
    /*======================*
     * 2) Calculate Eta obs *
     *======================*/
     
     Sequence.eta_obs = Calc_Eta_NSE(&Sequence);
    
                            
    //*========================================*
    //* 2) Iterate backwards to calculate ish  *
    //*========================================* 
    
    i=aa_count-1; //Not sure why we start at aa_count - 2, but this gives us the same eta variance
                  //calculated in msaum's semppr_2 code (originally aa_count-1 in msaum's code)
                  //8-14-13  whewgley - Changed back to aa_count - 1... the values calculated here deviated from
                  //msaums' semppr_2 code for some reason, but I changed this back to aa_count - 1 and now they are
                  //equivalent again. This seems to rely on whether the fasta read by semppr.data.load contains 
                  //a stop codon or not. 
                  
    aa_index_i = Sequence.aa_index[i];
    
    prod_e[i]=sq_prod_e[i]=prod_e_sq[i]=1; //by definition
    
    Y[i] = prod_e[i]*AA[aa_index_i].e_invc;
    varY[i] = AA[aa_index_i].e_invc2*prod_e_sq[i] - pow(AA[aa_index_i].e_invc*prod_e[i],2);
    H[i] = AA[aa_index_i].e_cpbinvc2*prod_e_sq[i]/(prod_e[i]*AA[aa_index_i].e_cpbinvc);
    
    i--;
    
    while(i>0){
        aa_index_ip1 = aa_index_i;
        aa_index_i=Sequence.aa_index[i];
        
        prod_e[i] = prod_e[i+1]*AA[aa_index_ip1].e_cpbinvc;
        prod_e_sq[i] = prod_e_sq[i+1]*AA[aa_index_ip1].e_cpb2invc2;
        sq_prod_e[i] = sq_prod_e[i+1]*AA[aa_index_ip1].e_cpbinvc*AA[aa_index_ip1].e_cpbinvc;
        
        Y[i] = prod_e[i]*AA[aa_index_i].e_invc;
        
        varY[i] = AA[aa_index_i].e_invc2*prod_e_sq[i] - pow(AA[aa_index_i].e_invc*prod_e[i],2);
        H[i] = AA[aa_index_i].e_cpbinvc2*prod_e_sq[i]/(prod_e[i]*AA[aa_index_i].e_cpbinvc);
        
        
        i--;
    }
    
    //Assume \sum_{j=n+1}^n(a_1+a_2j)Cov(Y_i,Y_j) = 0
    
    //*===================================================*
    //* 3) Iterate backwards again to calculate more ish  *
    //*===================================================*
    
    inner_loop_sum=0;
    var_eta=0;
    
    i = Sequence.aa_count - 1;
    
    a1pa2i = A1 + A2*(i);
      
    var_eta+=a1pa2i*a1pa2i*varY[i]+2*a1pa2i*Y[i]*inner_loop_sum;
    mean_eta+=a1pa2i*Y[i];
        
    inner_loop_sum+=a1pa2i*(H[i]-Y[i]);
    
    i--;
    
    double first_term=0,second_term=0,sum_first_term=0,sum_second_term=0;
    
    while(i>0){
        a1pa2i-=A2;
        first_term = a1pa2i*a1pa2i*varY[i];
        sum_first_term+=first_term;
        
        second_term = 2*a1pa2i*Y[i]*inner_loop_sum;
        sum_second_term+=second_term;
        
        var_eta+=a1pa2i*a1pa2i*varY[i]+2*a1pa2i*Y[i]*inner_loop_sum;
        mean_eta+=a1pa2i*Y[i];
        
        inner_loop_sum+=a1pa2i*(H[i]-Y[i]);
        
        i--;
    }
    
    var_eta = var_eta*B*B;
    mean_eta = B*mean_eta + A1+A2*(aa_count);
    
    //*==========================*
    //* 4) Find min and max eta  *
    //*==========================*
    
    eta_min_max(Sequence.aa_index,Sequence.aa_count,answer);
    
    //*========================================*
    //* 5) Fill in variables to pass back to R *
    //*========================================*
    
    *ETA_VAR = var_eta;
    *ETA_MEAN = mean_eta;
    *ETA_MIN = *(answer);
    *ETA_MAX = *(answer + 1);
    *ETA_OBS = Sequence.eta_obs;
                         
}

/*=======================================================================*
 *                      User Defined Functions                           *
 *=======================================================================*/
 
/*======================================================================================================*
 * FUNCTION: Process_R_input                                                                            *
 *      Purpose: This function is ONLY used in the sequence evolution                                   *
 *               code! It takes inputs from .C() R function (see                                        *
 *               http://users.stat.umn.edu/~geyer/rc/) and imports them                                 *
 *               to their corresponding variable names used in                                          *
 *               simulation_R-ext.c                                                                     *
 *      Inputs: CODON_INDEX - Codon Index Vector - vector containing integer values                     *
 *                respresenting each codon in the given sequence - Use codon index values               *
 *                found in preston/data/codon_summary.tsv                                               *
 *              ELONGATION_RATES - Codon Elongation Vector - vector containing codon elongation         *
 *                rates for each codon in order of codon indices                                        *
 *              MUTATION_RATES - Mutation Rate Vector - vector containing mutation rates                *
 *                for each codon in order of codon indices                                              *
 *              PHI - Phi value for the given codon sequence                                            *
 *              RESIDENCE_TIMES_SIM_GENOTYPE - Initialized vector to fill in the time each simulation   *
 *                step takes- must have length equal to number of evol steps                            *
 *              DELTA_ETA - Initialized vector to fill in (eta_sim - eta_obs) for each step             *
 *              UNSCALED_PR_SIM_GENOTYPES - Initialized vector to fill in \mu*exp{-y(eta_sim - eta_obs} *
 *                at each simulation step                                                               *
 *              UNSCALED_PR_OBS_GENOTYPE - Variable to hold the numerator of the likelihood calculation *
 *                for the original sequence                                                             *
 *              ETA_REF - Variable to hold the reference value for scaling eta to avoid                 *
 *                problems with numerical precision - This should be 0 if first                         *
 *                simulation for sequence, and for subsequent runs, this should                         *
 *                be a value that is passed back from the first simultion.                              *
 *                Now, the first simulation will pass back the first value of eta                       *
 *                as the reference value, so ETA_REF = ETA_OBS                                          *
 *              CC - Codon Counts Flag- this is unimportant to the current functionality                *
 *                and it will be removed NOTE *** REMEMBER TO REMOVE THIS                               *
 *              PO- Print Output Flag- this is unimportant to the current functionality                 *
 *                and it will be removed NOTE *** REMEMBER TO REMOVE THIS                               *
 *              POP - Effective Population Size- Parameter measuring genetic drift                      *
 *                DEFAULT = 1 (used to be 1.36*10^7 when dealing with phi)                              *
 *              A_1 - Ribosome initation cost in ATPs                                                   *
 *                DEFAULT a1=4                                                                          *
 *              A_2 - Ribosome elongation cost in ATPs, a2 in SEMPPR                                    *
 *                DEFAULT a2=4                                                                          *
 *              AT_BIAS - AT-Bias- Parameter related to mutation-- unused when including                *
 *                mutation file.                                                                        *
 *                DEFAULT = 0.5                                                                         *
 *              BENCH - Benchmark Flag- Flag that replaces rng numbers with predetermined               *
 *                numbers to set repeatable benchmarks                                                  *
 *                DEFAULT = 0                                                                           *
 *              BEE - Background Nonsense Error Rate B                                                  *
 *                DEFAULT = 0.0025                                                                      *
 *              GMT - Global Max Time- When the model is based on fixed time of evolution,              *
 *                this parameter limits the amount of time sequences are allowed to evolve.             *
 *                If the model is based on fixed number of evolution steps, this parameter              *
 *                does nothing.                                                                         *
 *                DEFAULT = -5                                                                          *
 *              IGNORE - When IGNORE is set equal to a positive integer n, nonsense errors can occur    *
 *                on the last n amino acids in a sequence and still produce a functional protein        *
 *                DEFAULT = 0                                                                           *
 *              RANDOM - Random Start Flag- Specifies whether the model should start with a random      *
 *                synonymous sequence or with the given sequence.                                       *
 *                DEFAULT = 0 (no random start)                                                         *
 *              GAMMA -Transition/Transversion Ratio- Parameter related to mutation-- unused when       *
 *                including mutation file.                                                              *
 *                DEFAULT = 1                                                                           *
 *              SCALE_FACTOR - Scaling Factor q- used to scale Effective Population Size Ne and         *
 *                Expression Rate phi.                                                                  *
 *                DEFAULT = 1 (used to be 4.18*10^-7 when dealing with phi)                             *
 *              MES - Minimum Evolution Steps- When the model is based on fixed number of evolution     *
 *                steps (substitutions), this paramter sets the number of evolution steps each          *
 *                sequence will be allowed.                                                             *
 *                DEFAULT = 10000                                                                       *
 *              BIS - Number of burn in steps per amino acid                                            *
 *                DEFAULT = 20                                                                          *
 *      Output: NONE
 *      Usage: refer to source/simulation_R-ext.c
 * =====================================================================================================*/

void Process_R_input(int *CODON_INDEX, double *ELONGATION_RATES, 
             double *MUTATION_RATES, int *AA_COUNT, double *PHI, 
             double *POP, double *A_1, double *A_2,
             double *AT_BIAS, double *BEE,
             int *IGNORE, double *GAMMA,
             double *SCALE_FACTOR, char ** AA_VEC,
             char ** CODON_VEC, int *N_CODON){
    int i, j, counter;
    //char tRNA_template[100]; //this is used to set up AA structs by reading in a blank trna file.
                         //the blank values are then filled in by values read from R


    // Read tRNA template file
    //TODO: 6/17/13 -- This information should be passed by R, not read in from a file here!
    //                   Information needed: Codon nts, aa codon corresponds to.
    //                   If ELONGATION_RATES and MUTATION_RATES order differs from AA, we're screwed!
    //sprintf(tRNA_template,"../preston/S.cerevisiae.tRNA.template.tsv");
    //n_aa = Read_tRNA_File(tRNA_template);
    
    //*=============================================*
    //*               STRUCT SETUP                  *
    //* Initialize AA structure and match codons to *
    //* amino acids                                 *
    //*=============================================*
    
    i=0,j=0;
    char currAA;
    int nCodons=0;
    char aa_list[22] = { 'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'Z', 'S', 'T', 'V', 'W', 'Y', 'X' };    //'X' represents the 'stop' amino acid, 'Z' represents the smaller set of serine codons: AGT and AGC.  The other ser codons are TCT, TCC, TCG, and TCA.
    
    //* Process info from R vecs *//
    
    
    //* Need to re-initialize num aa and num codons after every call *//
    n_aa = 0;
    for(i=0;i<22;i++){
		AA[i].num_codons = 0;
	}
    
    i = 0;
    while(i < *N_CODON){
        j=0;  
        currAA = AA_VEC[i][0];     
        AA[n_aa].aa = currAA;
        
        while(i<*N_CODON&&*(AA_VEC+i)[0]==currAA){
            strcpy(AA[n_aa].codon[j],CODON_VEC[nCodons++]);     //transfer codon string from R vector
            AA[n_aa].codon[j][3] = '\0';                        // must end with this character
            
            AA[n_aa].num_codons++;
            
            j++,i++;
        }
        
        n_aa++;
    }
    
    //* Check for errors *//

    if(n_aa >22) {
        Rprintf("Number of amino acids is greater than maximum value of 22\n");
        exit(EXIT_FAILURE);
    }
    
    for(i=0;i<n_aa;i++){
        j=0;
        
        if(AA[i].num_codons > 6){
            Rprintf("Number of codons for amino acid %d: %c exceeds maximum number of 6\n",i,AA[i].aa);
            exit(EXIT_FAILURE);
        }
        
        while(AA[i].aa!=aa_list[j]&&j<22){
            j++;
        }
        if(j == 22){
            Rprintf("Amino Acid %d: %c did not match any know amino acids\n",i,AA[i].aa);
        }
    }
    
    //*=============================================*
    //*             GENOME PARAMTERS                *
    //* These parameters apply to the entire genome *
    //*=============================================*
    
    counter = 0;
    for(i=0;i<n_aa;i++) {
        for(j=0;j<AA[i].num_codons;j++) {
            AA[i].elong_rate[j] = *(ELONGATION_RATES + counter);
            AA[i].incoming_mu[j] = *(MUTATION_RATES + counter);
            counter++;
        }
    }
    
    
    
    Ne = *POP;
    A1 = *A_1;
    A2 = *A_2;
    at_bias = *AT_BIAS;
    B = *BEE;
    gamma_ratio = *GAMMA;
    Q = *SCALE_FACTOR;
    ignore_aa = *IGNORE;
    
    //*==============================================*
    //*             SEQUENCE PARAMETERS              *
    //* These parameters only apply to this sequence *
    //*==============================================*
    
    Sequence.aa_count = *AA_COUNT;
    Sequence.phi_obs = *PHI;
    for(i = 0;i < *AA_COUNT;i++) {
        Sequence.codon_index[i] = *(CODON_INDEX + i);
    }

    /*
       Rprintf("POP=\t%f\n",*POP);
       Rprintf("A_1=\t%f\n",*A_1);
       Rprintf("A_2=\t%f\n",*A_2);
       Rprintf("AT_BIAS=\t%f\n",*AT_BIAS);
       //Rprintf("BENCH = %d\n",*BENCH);
       Rprintf("BEE=\t%f\n",*BEE);
       Rprintf("IGNORE=\t%d\n",*IGNORE);
       Rprintf("GAMMA=\t%f\n",*GAMMA);
       Rprintf("SCALE_FACTOR=\t%fe-7\n",*SCALE_FACTOR*10000000);
     */

}

void eta_min_max(int *aa_index,int aa_count, double * answer){
    // Purpose: calculate min and max eta, corresponding to optimal and pessimal sequences
    // Inputs: aa_index: memory address of first amino acid index 
    //         aa_count: number of amino acids in sequence
    // outputs: answer: memory address pointing to value of eta_min and eta_max
    //                  --note: *(answer) is eta min
    //                          *(answer + 1) is eta max
    
    int i,j;
    int pess_codon_index[2*aa_count]; //pessimal codon index (I allocated more than we really need)
    memset(pess_codon_index, 99, 2*aa_count*sizeof(int)); // Must explicitly initialize since aa_count is unknown to compiler
    int opt_codon_index[2*aa_count]; //optimal codon index
    memset(opt_codon_index, 99, 2*aa_count*sizeof(int));
    double t_sigma_vec_min[2*aa_count],t_sigma_vec_max[2*aa_count]; //temporary vec for elong pr's
    memset(t_sigma_vec_min, 1, sizeof(double));
    memset(t_sigma_vec_max, 1, sizeof(double));
    double t_xi,t_sigma_n; //more temp variables
    double t_min_elr,t_max_elr; //temporary variables for min and max elong rate
    int t_aa_index,t_codon_index; //temp var for aa_index
    double eta_min,eta_max;
    
    struct seq_struct optSeq;
    struct seq_struct pessSeq;
    
    //*============================================*
    //* 1) Find optimal and pessimal codon indices *
    //*============================================*
    
    for(i=0;i<aa_count;i++){
        t_aa_index = *(aa_index+i);
        t_min_elr = AA[t_aa_index].min_elong_rate;
        t_max_elr = AA[t_aa_index].max_elong_rate;
        
        //* Search through codons for t_min_elr and t_max_elr
        for(j=0;j<AA[t_aa_index].num_codons;j++){           
            t_codon_index = AA[t_aa_index].codon_index[j];
            
            //* If this codon has AA's min elong rate, store in pessimal
            //* sequence
            if(Codon[t_codon_index].elong_rate == t_min_elr){
                pessSeq.codon_index[i] = t_codon_index;
            }
            
            //* If this codon has max elong rate, store in optimal
            if(Codon[t_codon_index].elong_rate == t_max_elr){
                optSeq.codon_index[i] = t_codon_index;
            }
        }
        
        //* Check to see if this aa position was assigned an optimal
        //* and pessimal codon index
        
        if(pessSeq.codon_index[i]==25443||optSeq.codon_index[i]==25443){ //25443 is the integer equivalent of the character 99
            Rprintf("Error in eta_min_max(): AA position %d was not given an optimal or pessimal codon index",i);
        }
    }
    
    
    //*==========================*
    //* 2) Calculate minimum eta *
    //*==========================*
    
    eta_min = Calc_Eta_NSE(&optSeq);
        
    //*==========================*
    //* 3) Calculate maximum eta *
    //*==========================*
        
    eta_max = Calc_Eta_NSE(&pessSeq);
        
    //*===================================*
    //* Store eta min and max in answer[] *
    //*===================================*
    
    answer[0] = eta_min;
    answer[1] = eta_max;
        
}

void Codon_Counts(int * codon_index_vec ,int * codon_cts_vec ,int aa_count){
    int j;
    int codon_indx;
    
    //* Initialize
    for(j=0;j<61;j++){
        *(codon_cts_vec+j)=0;
    }
    
    //* Add up codon counts
    for(j=0;j<aa_count;j++){
        codon_indx = *(codon_index_vec + j);
        codon_cts_vec[codon_indx]++;
    }
}

void Rprint_codon_sequence( int *codon_index ){
    int i;
    Rprintf("Codon Sequence Start:\n");
    for(i=0;i<Sequence.aa_count;i++){
        Rprintf("%d\t",*(codon_index + i));
    }
    Rprintf("\n\n");
}



void bin_eta(double eta, double * bins, double * bin_lims, double value, int num_bins){
    //Purpose: Place eta into bin for histogram. This function will look for the bin interval
    //         (determined by bin_lims) that eta fits into and add the desired number 
    //         (value) to that bin.
    //Inputs: eta: eta value of interest
    //        bin_lims: interval limits for bins. must be in increasing order
    //        value: value added to bin 
    //        num_bins: number of bins, equal to length of bin_lims - 1
    //Outputs: None
    
    int i=0;
    
    
    //*=============================*
    //* Search for eta in bin_lims  *
    //*=============================*
    
    if(eta < *(bin_lims)){
        error("Observed Eta: %f less than min bin limit: %f\n",eta,*(bin_lims));
        
    }
    if(eta > *(bin_lims + num_bins)){
        error("Observed Eta: %f greater than max bin limit %f\n",eta, *(bin_lims+num_bins));
    }
    
      //* Eta should be between two values of bin_lims
    while((!(eta >= *(bin_lims + i-1) && eta < *(bin_lims + i)))&&i<=num_bins) i++;
    
    i--; //counteract i++
    
    //*===================*
    //* Add value to bin  *
    //*===================*
    //Rprintf("eta_bin_low: %f eta: %f eta_bin_high: %f\n",*(bin_lims + i),eta,*(bin_lims + i+1));
    *(bins + i) += value; 
        
}

void intervals(double from, double to, int length, double * int_lims){
    //Purpose: create array sequence of length 'length' from 'from' to 'to' 
    
    int i;
    double interval;
    
    interval = (to - from)/(length-1);
    
    *(int_lims) = from;
    
    for(i=1;i<length;i++){
        *(int_lims + i) = from + interval*i;
    }
    
    //check if last entry equal to 'to'
    if(*(int_lims + length-1)!=to) error("Sequence function did not work\nLast entry: %f\tto: %f",*(int_lims+length-1),to);
    
    
}


 void flip_codons(struct seq_struct *currSeq, struct seq_struct *propSeq, double *randNums, int num_flips){
	 int i,j;
	 int * currCodonIndex = currSeq->codon_index;
	 int currCodon,currAA;
	 int seqLen = currSeq->aa_count; //length of current sequence
	 int flip_count = 0;
	 int newCodon;
	 
	 while(flip_count < num_flips){
		 
		 i = j = 0;
		 //1) Find the codon index that needs flipped
		 while(i < *(randNums + 2*flip_count)*seqLen&&i<seqLen) i++;
		 i--;
		 
		 currCodon = *(currCodonIndex+i);
		 currAA = Codon[currCodon].aa_index;
		 
		 //2) Find the codon it will be flipped to
		 while(j < *(randNums + 2*flip_count+1)*AA[currAA].num_codons) j++;
		 j--;
		 
		 newCodon = AA[currAA].codon_index[j];
		 
		 //3) Change codon in propSeq
		 propSeq->codon_index[i] = newCodon;
		 
		 
	
		 flip_count++;
	} 
	 
 }
