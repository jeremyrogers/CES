/*==========================*
 * Preprocessor directives  *
 *==========================*/
#ifndef _CES_struct
#define _CES_struct
 
#define MAX_DDIM 30000
#define MAX_AA 6500
#define MAX_NTS 19500
#define MAX_LOCI 6000		//set max loci to read in and simulate
#define MAX_TIME -5		//Max number of evolutionary steps to simulate,
		    // if set to a negative value it will calculate the appropriate max time to have X substitutions/nt.
		    //no longer used.  Defined via coommand line
#define NE 1.36e7		//effective population size.
		  //#define Q 1.519e-5 //fitness scaling coefficient
#define MU 1E-9			//per generation mutation rate


/*=======================================================*
 *  STRUCTURE: amino_acid								 *
 * 														 *
 *	   Defining structure for a single amino acid        *
 *	   keep track of E(1/c), E(1/c^2), E((c+b)/c)... 	 *
 *	   which are expectations wrt the translation rates  *
 *	   for the aa's set of codons 						 *
 * ======================================================*/
 
struct amino_acid {
      char codon[6][4];
      char aa;
      int num_codons;		//was cc
      int codon_index[6];	//used when generating random sequences

      //expected values for elongation related rates
      double e_invc;				// 1/c
      double e_cpbinvc;				//(c+b)/c
      double e_invc2;				//1/c^2
      double e_cpb2invc2;			//(c+b)^2/c^2
      double e_cpbinvc2;			//(c+b)/c^2 --used in calc cov b/w codons
      double elong_rate[6];			//elongation rate of codons for AA, was tr_rate
      double elong_pr[6];			//Pr of successful elongation c/(c+b)
      double fail_pr[6];			//Pr of unsuccessful elongation b/(c+b)
      double max_elong_rate;		//max elongation rate for the amino acid
      double max_elong_pr;			//max elongation rate for the amino acid
      double min_elong_rate;		//min elongation rate for the amino acid
      double min_elong_pr;			//min Pr of successful elongation for the amino acid
      double neutral_obs_pr[6];		//Pr of observing this codon based solely on its nt sequence and any AT bias.
      double neutral_obs_pr_cum[6];	//Cumulative probability of observing this codon. Useful for when choosing a random codon when there is AT bias
      double mu[6][6];				//hold mutation values
      double incoming_mu[6]; 					  		//sum of all mutation rates to codon
};



/*=======================================================*
 * STRUCTURE: codon_struct								 *
 * 														 *
 * 		Defining structure for a single codon. 			 *
 * 		Keeps track of individual information as well as *
 * 		information relative to one-step neighbors.		 *
 * ======================================================*/

struct codon_struct {
      char codon[4];
      char aa;			//aa letter

      int synonym_index[5];	//Alternative synonyms for the same aa listed by their codon index #
      int one_step_synonym_index[5];	//as above but subset that is one step mutant from codon.

      int aa_index;

      int num_synonym;
      int num_one_step_synonym;
      


      double elong_pr_ratio[6];						//(c_i/(c_i+b))/(c_j/(c_j+b)) for one step mutant from codon c_i.
      double delta_b_over_c[6];						//(b/c_i - b/c_j) for one step mutant from codon c_i.
      double relative_mu; 							// relative mutation rate to this codon from all one step synonyms
      double one_step_relative_mutation_rate[5];	//mutation rate to one step neighbors of current codon
													//as above but subset that is one step mutant from codon.
	  double syn_relative_mutation_rate[5];    //mutation rate from this codon to each neighbor (all neighbors,
	                                                //not just one-step neighbors)
	  double incoming_mu;
	 

      double b_over_c;		//b/c_i
      double elong_rate;	//elongation rate of codons for AA, was tr_rate
      double elong_pr;		//Pr of successful elongation c/(c+b)
      double fail_pr;		//Pr of unsuccessful elongation b/(c+b)
};



/*==========================================================*
 * STRUCTURE: seq_struct									*
 * 															*
 * 		Defining structure for one genetic sequence.		*
 * 		Holds information relevant to the whole sequence	*
 * =========================================================*/

struct seq_struct{
      char id[26];		//Gene ID or ORF name
      char codon_seq[MAX_AA][4];	//codon sequence
									
      int codon_index[MAX_AA];	//codon index code
      int aa_index[MAX_AA];		//aa index code

      int aa_count;						//aa count
      double sigma_vec[MAX_AA];			//vector of \sigma(i) values.  
										//Note that sigma_0 repesents the start codon so its value is 1. 
      double sigma_ratio_vec[MAX_AA];	//vector of \sigma(i)/\sigma(n) values
      double max_time;
      int ARR;
      double pA;		//alpha parameter for f(eta)~beta distribution
						//used by gsl
      
      
	  /* Dont need any of these
	  
      double pB;		//beta+alpha parameter for f(eta)~beta distribution
      double pD;		//eta_diff = eta_max - eta_min
      double pN;		//eta_obs-eta_min
      double pdenom;		//Scaling term for pdf
      double MEAN;		//arithmetic mean expression ?
      double z_max;
      double z_mode;		//mode of z based on eta~beta dist
      double z_prcntl[202];
      double z_amean;
      double z_gmean;
      double z_var;
      double z_sd;
      double z_mode_gamma;	//mode based on assuming eta ~ gamma
      */
      double sigma_obs;			//obs sigma value
      double xi_obs;			//obs xi value
      double eta_obs;			//obs eta value
      double mu_obs;            //mutation coefficient for observed sequence = prod_i \mu_i^(x_i(\cvec))
                                //  where mu_i is the sum of relative mutation rates to codon i and x_i
                                //  is the number of codon i in the sequence cvec  
      
      /*
      double eta_min;		//min eta value
      double eta_max;		//max eta value
      double eta_mean;		//mean of eta distn
      double eta_var;		//var of eta distn
      double eta_sd;		//sd of eta distn

      double alpha;		//alpha parameter for f(eta)~beta distribution
      double beta;		//beta parameter for f(eta)~beta distribution
      double gamma_alpha;	//alpha parameter for f(eta)~gamma distribution
      double gamma_beta;	//beta parameter for f(eta)~gamma distribution
	  */
	  int codon_cts[64]; //running total of number of codons
	  
      // following added for evolution simulations
      double eta_initial;		//eta value at start of simulation
      double phi_obs;			//assumed production rate
      double z_obs;				//assumed production rate, scaled by Ne and q
      double delta_eta_mean;	//mean of distn of delta_eta values for wt allele
      double delta_eta_var;		//var of distn of delta_eta values for wt allele

      double evol_delta_eta_mean;	//mean diff b/w wt_eta and mut_eta;
      double evol_delta_eta_var;	//var b/w wt_eta and mut_eta;
      double evol_time;				//var b/w wt_eta and mut_eta;

      int evol_steps;		//number of times wt allele is replaced by a mutant
      int num_synon_mut;	//num of mut that were synonymous

};

#endif /* _CES_struct */
