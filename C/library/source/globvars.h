/*==============================*
 *  Declaring Global Variables  *
 *==============================*/
 
extern double A1;				//Cost of initiation of translation in units of ATP. Default val set in command line to  A1 = 2
extern double A2;				//Cost each elongation step in ATPs, default set in command line to A2=4
extern double B;				//Nonsense error rate. units of 1/sec. Default value set in commandline function to   B = 0.00515
extern double pi;				//
extern double mu;				//per nt mutation rate.
extern double min_z_max;	    //minimum value to start integration routines with;
extern double max_z_max;	    //set an upper bound for z
extern double z_max_factor;	    //amount to multiply z_max by when searching for value
extern double relerr;		    //relative error for integration routines
extern double abserr;		    //absolute error for integration routines.  aka tolerance
extern double at_bias;		    //parameter for adjusting nt composition due to biased mutation and/or gene conversion.
							    //According to table 6.1 in Lynch (2007) on p. 125, the observed AT bias is 0.62
extern double mu_bias;			//parameter calculated based on at_bias.  mu_bias = at_bias/(1-at_bias)
extern double gamma_ratio;		//ratio of transition mutation rate to transversion mutation rate 
							    //This is \alpha/\beta in the Tamura and Nei (1993) where we assume \alpha_1 = \alpha_2
							    //Calculations using data from Lynch give value of approximately 1.22 for S.c. 
extern double mutation_matrix[4][4];	    //mutation matrix order is ATCG.  See mutation.matrices.nb for details.

extern double global_max_time;	//Max number of evolutionary steps to simulate, if set to a negative value it will set the max time so that, under neutrality, there are on average X substitutions/nt.
extern double Ne;//NE;
extern double Q;//4.19e-7;		            // Fitness scaling coefficient

extern int ignore_aa;				        //decrease true aa_count by this amount.
extern int random_start;		            //indicate whether or not to start with a random sequence (1) or the seq read in (0);
extern int n_aa,				            //Number of aa in a sequence
    n_seq,				                    //set when reading fasta file
    pout,					                //print out configuration flag.  Stores -W arguement passed on command line
    pconf, gconf,
    runs, analytic, max_aa, MLE;
extern int print_evol_fasta;		        //indicate whether or not to print seq evolution to fasta file
                                            //*cannot* be set via command line
extern int print_delta_eta_vals;	        //indicate whether or not to print list of delta eta vals for each step
								            //can be set to 1 via command line
extern int print_eta_trace;		            //indicate whether or not to print the eta values and residence times for the wild type genotypes
								            //can be set to 1 via command 
extern int print_debug;								            

extern int benchmark;			        	//Uses predetermined values to generate random sequence and calculate mut_pr and time_step
								            //instead of random generated 
extern int codon_counts;

extern int simulation_steps;
extern int burn_in_steps;

extern int non_random_time_step;    //If this is equal to 0, time step will be pulled from a random exponential dist, if this 
								    //is equal to 1, time step will be calculated based on leaving rate	**DEFAULT 1**
extern int mutation;

extern int R_or_C;     //This flag is used to determine whether the CES library functions should use printf or Rprintf
				       //R_or_C = 0 -> Calling from R, use Rprintf
