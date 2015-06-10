/*===========================*
 *  Declaring Header Files   *
 *===========================*/

#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>		//for mkdir
#include <sys/types.h>		//for types w/in mkdir
#include <math.h>
#include <time.h>
#include <sys/time.h>		//gettimeofday()
#include <string.h>
#include <omp.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_min.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_rng.h>	//uniform rng
#include <gsl/gsl_randist.h>	//rng from distns

/*==================================*
 * Include CES library header files *
 *==================================*/

#include "CES.h"

/*=======================*
 * Turn off error codes  *
 *=======================*/

int gsl_errors_off = 0;
int print_debug = 0;

/*==============================*
 *  Declaring Global Variables  *
 *==============================*/

double A1;					//Cost of initiation of translation in units of ATP. Default val set in command line to  A1 = 2
double A2;					//Cost each elongation step in ATPs, default set in command line to A2=4
double B;					// Nonsense error rate. units of 1/sec. Default value set in commandline function to   B = 0.00515
double pi = M_PI;				//
double mu = MU;				//per nt mutation rate.
double min_z_max = 1e-10;	//minimum value to start integration routines with;
double max_z_max = 1e100;	//set an upper bound for z
double z_max_factor = 10;	//amount to multiply z_max by when searching for value
double relerr = 1e-6;		//relative error for integration routines
double abserr = 1e-20;		//absolute error for integration routines.  aka tolerance
double at_bias = 0.5;		//parameter for adjusting nt composition due to biased mutation and/or gene conversion.
//According to table 6.1 in Lynch (2007) on p. 125, the observed AT bias is 0.62
double mu_bias = 1;			//parameter calculated based on at_bias.  mu_bias = at_bias/(1-at_bias)
double gamma_ratio = 1;		//ratio of transition mutation rate to transversion mutation rate 
//This is \alpha/\beta in the Tamura and Nei (1993) where we assume \alpha_1 = \alpha_2
//Calculations using data from Lynch give value of approximately 1.22 for S.c. 
double mutation_matrix[4][4];	//mutation matrix order is ATCG.  See mutation.matrices.nb for details.

double global_max_time = MAX_TIME;	//Max number of evolutionary steps to simulate, if set to a negative value it will set the max time so that, under neutrality, there are on average X substitutions/nt.
double Ne = 1;//NE;
double Q = 1;//4.19e-7;		// Fitness scaling coefficient

int ignore_aa;				//decrease true aa_count by this amount.
int random_start = 0;		//indicate whether or not to start with a random sequence (1) or the seq read in (0);
int n_aa = 0,				//Number of aa in a sequence
		n_seq,				    //set when reading fasta file
		pout,					//print out configuration flag.  Stores -W arguement passed on command line
		pconf, gconf,
		//  calcamean=0,        //indicates whether to calculate arithmetic mean.  Will be reset based on arguments passed
		//  calcgmean=0,
		//  calcvar=0,
		//  calcmode = 1,       //indicates which method to use to estimate the posterior mode. 1:gamma based 2: numerical search
		runs, analytic, max_aa, MLE;
		int print_evol_fasta = 0;		//indicate whether or not to print seq evolution to fasta file
		//*cannot* be set via command line
		int print_delta_eta_vals = 0;	//indicate whether or not to print list of delta eta vals for each step
		//can be set to 1 via command line
		int print_eta_trace = 0;		//indicate whether or not to print the eta values and residence times for the wild type genotypes
		//can be set to 1 via command line

		int benchmark = 0;				//Uses predetermined values to generate random sequence and calculate mut_pr and time_step
		//instead of random generated 
		int codon_counts = 0;

		/*=========================================*
		 *     Code Specific Global Variables      *
		 *  These variables are different in the   *
		 *  data simulation and sequence evolution *
		 *  code.                                  *
		 *=========================================*/

		int non_random_time_step = 0;   //If this is equal to 0, time step will be pulled from a random exponential dist, if this 
		//is equal to 1, time step will be calculated based on leaving rate **DEFAULT 0**
		int mutation = 0;
		int R_or_C = 1; //This flag is used to determine whether the CES library functions should use printf or Rprintf
		//R_or_C = 1 -> C only, use printf

		/*=====================*
		 * default input files *
		 * ====================*/

		char def_tRNA[150] = "tRNA_link";
		char fasta_file[150] = "fasta_link";
		char out_prefix[150], trna[150], mutation_file[150] = "NONE";
		char out_folder[150];		//derived from out_prefix
		char command_line[1000];	//save command line for printing
		char exec_time[360];		//save time command executed at.



		time_t start;
		struct tm *timeinfo;
		struct timeval time_val, comp_start, comp_stop, print_stop, time_diff;

		/*================*
		 * Define Structs *
		 *================*/

		struct amino_acid AA[22];
		struct codon_struct Codon[64];
		struct seq_struct Seq[6000];

		/*======================================================*
		 * FUNCTION: Read_Phi_File								*
		 * 		This function reads a file holding protein		*
		 * 		production level values.						*
		 * 		See readme for details concerning file formats	*
		 * =====================================================*/

		void
		Read_Phi_File (char *filename, int num_seq) {
			FILE *file_handle;
			char curr_char;
			char seq_id[26];
			double phi;
			int i = 0, j = 0;
			int fscanf_check;

			//Initialize phi values
			for (i = 0; i < num_seq; i++) {
				Seq[i].phi_obs = 0;
			}

			file_handle = fopen (filename, "r");
			if (!file_handle) {
				printf ("\nPhi File: %s Doesn't Exist\n", filename);
				wrong ();
			}

			curr_char = fgetc (file_handle);

			//If header, read untill new line.
			if (curr_char == '"') {
				while (curr_char != '\n') {
					curr_char = fgetc (file_handle);
				}
				curr_char = fgetc (file_handle);
			}
			//This while loop reads til end of file
			while (curr_char != EOF) {
				i = 0;
				j = 0;
				//This while loop reads sequence ID     
				while (curr_char != ',' && curr_char != '\t') {
					seq_id[i] = curr_char;
					curr_char = fgetc (file_handle);
					i++;
				}
				seq_id[i] = '\0';
				//Read in Phi value
				fscanf_check = fscanf (file_handle, "%lf", &phi);
				if (fscanf_check < 1) {
					fprintf(stderr, "Error reading phi value\n");
				}
				//Find Sequence with seq_id
				while (strcmp (seq_id, Seq[j].id) != 0 && j < num_seq)
					j++;
				//Store Phi in Seq[j]
				if (j != num_seq)
					Seq[j].phi_obs = phi;


				//Read til \n or EOF
				while (curr_char != '\n' && curr_char != EOF)
					curr_char = fgetc (file_handle);
				//Get next line if not EOF
				if (curr_char == '\n')
					curr_char = fgetc (file_handle);
			}

		}

/*======================================================*
 * FUNCTION: Read_Mutation_File_1						*
 * 		This function opens and reads an incoming sum	*
 * 		relative mutation rate file. To be done after	*
 * 		reading in tRNA file so we can check num codons	*
 * 		is correct. 									*
 * =====================================================*/

	void
Read_Mutation_File_1 (char *filename)
{
	int i, j, k;
	FILE *file_handle;
	char curr_char;
	char curr_codon[4];
	double mu_rate;
	int num_codons;
	int max_aa = 22;		//20 AA but serine may get split and we may have stop codons
	//int codon_index;
	int test;
	int fscanf_check;

	j = 0;

	file_handle = fopen (filename, "r");
	if (!file_handle) {
		printf ("\nMutation File: %s Doesn't Exist\n", filename);
		wrong ();
	}

	curr_char = fgetc (file_handle); //get aa index 


	while (curr_char != EOF) {
		i = 0;		//set AA index to 0

		while (curr_char != AA[i].aa && i < max_aa) {  //flip through AA until you get a match
			i++;
		}

		if (i == max_aa) {
			if (curr_char == '"') {
				while (curr_char != '\n'){
					curr_char = fgetc (file_handle);
				}
				//Get AA index for next line
				curr_char = fgetc (file_handle);
			}else {
				printf ("\nAA index '%c' in file did not match any known AA. Exiting...\n",curr_char);
				exit (1);
			}
		} else {
			curr_char = fgetc (file_handle);	//get next character. Should be a \t

			/* Read in codon */

			for (j = 0; j < 3; j++)	//load codon sequence and put in codon_index codon index
			{
				curr_codon[j] = fgetc (file_handle);
			}
			curr_codon[j] = '\0';


			curr_char = fgetc (file_handle);	//get next char. Should be a \t

			/* Read in mutation rate */

			fscanf_check = fscanf (file_handle, "%lf", &mu_rate);
			if (fscanf_check != 1) {
				fprintf(stderr, "Error reading mutation file\n");
			}

			/* Increment codon count for the aa */

			num_codons++;


			/*Read until end of line or EOF */

			while ((curr_char != '\n') && (curr_char != EOF))
				curr_char = fgetc (file_handle);

			/* Cycle through codons to find match */

			j = 0;
			test = 0;
			while ((test == 0) && (j < AA[i].num_codons)) {
				if (!strcmp (AA[i].codon[j], curr_codon)) {
					AA[i].incoming_mu[j] = mu_rate;	
					test = 1;	//codons match.  Exit loop
				} else {
					j++;
				}
			}
			if (j >= AA[i].num_codons) {
				printf ("\nCurrent codon %s did not match any codons for aa %d. Exiting...\n", curr_codon, AA[i].aa);
				exit (1);
			}

			/* Get aa index for next codon if not at EOF */

			if (curr_char != EOF)
				curr_char = fgetc (file_handle);

		}
	}

	/* ================================================================*
	 * Since we cannot find individual codon mutation rates, we assume *
	 * all incoming mutation rates from synonymous codons are equal.   *
	 * That is, for example, given a three codon case c_i,c_j,c_k each *
	 * with mu_i=(Sum of synonymous codon mutation rates  to codon i)= *
	 * mu(j to i)+mu(k to i), we say mu(j to i)=mu(k to i)=mu_i/2.     *
	 * Since stationary state depends only on mu_i, choosing arbitrary *
	 * values for synonymous codon mutation rates should conserve      *
	 * equilibrium. <whewgley> 7-26-12                                 *
	 * ================================================================*/

	for (i = 0; i < max_aa; i++) {
		for (j = 0; j < AA[i].num_codons; j++) {
			//mikeg: I don't think this should be divided by num_codons, but instead num_neighbors (one_step_neighbors)
			//       non-accessible synonymous codons should have mutation rates of zero.
			if(AA[i].num_codons > 1) mu_rate = AA[i].incoming_mu[j] / (AA[i].num_codons - 1);
			for (k = 0; k < AA[i].num_codons; k++) {
				if (j == k)
					AA[i].mu[k][j] = 0;
				else
					AA[i].mu[k][j] = mu_rate;
			}
		}
	}
	fclose (file_handle);
}


/*======================================================*
 * FUNCTION: Read_Mutation_File_2						*
 * 		This function opens and reads an individual		*
 * 		relative mutation rate file. To be done after	*
 * 		reading in tRNA file							*
 * =====================================================*/

void
Read_Mutation_File_2 (char *filename) {
	char curr_char = { '\0' }, col_codons[6][4] = {
		"\0\0\0\0"}, curr_AA, row_codon[4] = {
			'\0'};
	double mut_rate_matrix[6][6] = {{0}};
	double incoming_sum_mu;
	int i = 0, j = 0, num_codons, indx;
	int max_aa = 22;
	int codon_index;
	int fscanf_check;

	FILE *file_handle;

	file_handle = fopen (filename, "r");

	if (!file_handle) {
		printf ("\nMutation File: %s Doesn't Exist\n", filename);
		wrong ();
	}

	curr_char = fgetc (file_handle);

	/*Check for header: ignore for now WHEWGLEY */
	if (curr_char == '"') {
		while (curr_char != '\n')
			curr_char = fgetc (file_handle);
	} else if (curr_char == '>') {
		fseek (file_handle, 0, SEEK_SET);
	} else {
		printf ("Error in mutation file format 1.\n");
		exit (-1);
	}

	curr_char = fgetc (file_handle);

	while (curr_char != EOF) {



		i = 0;

		if (curr_char != '>') {
			printf ("Error in mutation file format 2.\n");
			exit (-1);
		}

		curr_char = curr_AA = fgetc (file_handle);	//should be first AA 
		while (AA[i].aa != curr_AA) {
			i++;
			if (i == max_aa) {
				printf ("Did not recognize Amino Acid '%c' in mutation file\n", curr_AA);
				exit (-1);
			}
		}

		indx = i;

		num_codons = AA[indx].num_codons;

		for (i = 0; i < num_codons; i++) {
			fscanf_check = fscanf (file_handle, "\t %s", col_codons[i]);
			if (fscanf_check < 1) {
				fprintf(stderr, "Error reading mutation file\n");
			}
		}

		curr_char = fgetc (file_handle);
		while (curr_char != '\n')
			curr_char = fgetc (file_handle);

		for (i = 0; i < num_codons; i++) {
			fscanf_check = fscanf (file_handle, "%s", row_codon);
			if (fscanf_check < 1) {
				fprintf(stderr, "Error: wrong mutation file format\n");
			}

			if (strcmp (row_codon, col_codons[i]) != 0) {
				printf ("Wrong mutation file format\n");
				exit (-1);
			}

			for (j = 0; j < num_codons; j++) {
				fscanf_check = fscanf (file_handle, "\t%lf", &mut_rate_matrix[i][j]);
				if (fscanf_check < 1) {
					fprintf(stderr, "Error: wrong mutation file format\n");
				}
			}

			curr_char = fgetc (file_handle);
			while (curr_char != '\n' && curr_char != EOF)
				curr_char = fgetc (file_handle);


		}

		for (i = 0; i < num_codons; i++) {
			for (j = 0; j < num_codons; j++) {
				AA[indx].mu[i][j] = mut_rate_matrix[i][j];
			}
		}

		incoming_sum_mu = 0; //re-initialize 

		for(j = 0; j < num_codons; j++) {
			for (i = 0; i < num_codons; i++) {
				incoming_sum_mu += AA[indx].mu[i][j];
			}
			AA[indx].incoming_mu[j] = incoming_sum_mu;
			codon_index = AA[indx].codon_index[j];
			Codon[codon_index].incoming_mu = incoming_sum_mu;
		}
		curr_char = fgetc (file_handle);
	}


	return;
}

/*======================================================*
 * FUNCTION: Read_FASTA_File							*
 * 		Reads fasta file into Seq structures.			*
 * 		Codons stored in Seq.codon_seq[] array.
 * =====================================================*/

int
Read_FASTA_File (char *filename, int *phi_flag) {

	int i, j, k;
	char curr_char;		//, codon[4];
	char curr_line[82];
	char *char_ptr;
	char *str_test;
	char nt_seq[MAX_NTS];
	int end_of_seq;
	int nt_count, aa_count;
	FILE *fh;

	fh = fopen (filename, "r");
	if (!fh) {
		printf ("\nFASTA/Sequence File Doesn't Exist\n");
		wrong ();
	}
	//codon[3]='\0';

	//  curr_char='';
	str_test = curr_line;

	//Read until you get to a new line with a >
	do {
		str_test = fgets (curr_line, 82, fh);
	}
	while (curr_line[0] != '>' && str_test != NULL);

	i = 0;

	while (str_test != NULL && i < MAX_LOCI) {

		//line should be start of new sequence
		//read in seq id
		k = 0;
		curr_char = curr_line[0];
		while ((curr_char != ' ') && (curr_char != '\t')
				&& (curr_char != '\n')
				&& (k < 26)) {
			curr_char = curr_line[k + 1];	//offset by 1 b/c of >
			Seq[i].id[k] = curr_char;
			k++;
		}
		if (k == 25) {
			printf ("Seq ID greater than 25 characters. Exiting\n");
			exit (1);
		}
		//replace new line with null character
		Seq[i].id[--k] = '\0';

		//Read in phi value
		//locate start of text  "phi = "
		char_ptr = strstr (curr_line, "phi=");	//removed space b/w phi and = sign 

		if (char_ptr == NULL)
			(*phi_flag) = 0;	/*Set flag if phi values are not included in fasta.
													Program will exit if phi_flag==0&& no -P option. */
		else
			(*phi_flag) = 1;

		if (*phi_flag == 1) {
			char_ptr += strlen ("phi =");	//go to end of match string
			//note strtod will skip any whitespace preceeding the number
			Seq[i].phi_obs = strtod (char_ptr, NULL);	//first argument is ptr to str location of double
			// arguementsecond is a pointer to the place afterwards
		}

		nt_count = 0;
		nt_seq[0] = '\0';


		end_of_seq = 0;


		//get next line
		str_test = fgets (curr_line, 82, fh);

		while (str_test != NULL && end_of_seq == 0) {
			//test to make sure string isn't too long
			if (strlen (curr_line) > 81) {
				printf ("ERROR: FASTA file column is >80 characters. Exiting");
				exit (1);
			}
			//read through until you hit a '>'
			curr_char = curr_line[0];

			switch (curr_char) {
				case 'A':
				case 'T':
				case 'G':
				case 'C':
				case 'a':
				case 't':
				case 'g':
				case 'c':
					k = strlen (curr_line);
					if (curr_line[k - 1] == '\n'
							|| curr_line[k - 1] == EOF)
						k--;
					nt_count += k;
					if (nt_count >= MAX_NTS) {
						printf ("Seq[%d].id = %s has %d nts.  MAX_NTS is %d. Exiting", i, Seq[i].id, nt_count, MAX_NTS);
						exit (1);
					}
					strncat (nt_seq, curr_line, k);
					//get next line
					str_test = fgets (curr_line, 82, fh);
					break;
				case '>':
					end_of_seq = 1;
					//don't get next line
					break;
				case '\n':
					//get next line
					str_test = fgets (curr_line, 82, fh);
					break;
				case EOF:	//str_test!=NULL should prevent reaching this point
					end_of_seq = 1;
					break;
				default:
					printf ("Read in unexpected character (%c) at beginning of line. Exiting.", curr_char);
					exit (1);
					break;
			}
		}

		//process nt seq
		if (nt_count % 3 != 0) {
			printf ("Seq[%d].id = %s had %d nts, which is not a multiple of 3. Exiting", i, Seq[i].id, nt_count);
			exit (1);
		} else {		//everything seems okay
			aa_count = (nt_count / 3) - 1;	//subtract 1 b/c of stop codon
			Seq[i].aa_count = aa_count;

			k = 0;

			for (j = 0; j <= aa_count; j++) {	//include stop codon
				Seq[i].codon_seq[j][0] = nt_seq[k++];
				Seq[i].codon_seq[j][1] = nt_seq[k++];
				Seq[i].codon_seq[j][2] = nt_seq[k++];
				Seq[i].codon_seq[j][3] = '\0';
			}

		}

		i++;		//increment seq id

	}
	return i;
}


/*======================================================*
 * FUNCTION: Convert_Codon_Seq_to_Codon_Index			*
 * 		Takes the codon_seq array from a Seq and 		*
 * 		assigns each codon a unique index (1-64) and 	*
 * 		stores it in a new vector index_vec.			*
 * =====================================================*/

void
Convert_Codon_Seq_to_Codon_Index (char codon_seq[][4],
		int *index_vec,
		int aa_count) {
	int i, j, match;
	//  char codon[4];

	for (i = 0; i < aa_count; i++) {

		//    strcpy(codon, codon_seq[i]);

		j = 0;
		do {
			match = strcmp (Codon[j++].codon, codon_seq[i]);
		}
		while (match != 0 && j < 64);

		if (j == 64) {
			printf ("ERROR: Can't match %s to a codon. Exiting",
					codon_seq[i]);
			exit (1);
		} else {
			index_vec[i] = --j;	//counter act j++ in do loop
		}

	}
}



/*======================================================*
 * FUNCTION: Read_Commandline_Args						*
 * 		Takes in options and reads input files			*
 * =====================================================*/

void
Read_Commandline_Args (int argc, char *argv[]) {
	char *dp;
	struct stat sb;		//for checking folder status
	int i, j, test;
	int phi_flag;		//sets a flag if fasta does not have included phi
	int ncheck = 19;		//should be equal to size of check[]
	int check[19] = { 0 };	//,0,0,0,0,0,0,0,0,0,0,0,0,0}; //to add more command line arguments increase this by 0
//trna, fasta_file, out

if (print_debug) {
	printf ("reading arguments\n");
	fflush (stdout);
}
for (i = 1; i < argc; i++) {
	if (argv[i][0] == '-') {
		switch (argv[i][1]) {
			case 'T':
			case 't':
				if ((argv[i][2] != '\0') || (i == argc - 1)) {
					printf ("\ntRNA file name not specified or Incorrect usage\n");
					wrong ();
				} else {
					//n_aa=Read_Files(argv[++i],1);
					n_aa = Read_tRNA_File (argv[++i]);
					strcpy (trna, argv[i]);
					check[0] = 1;
					break;
				}
			case 'F':
			case 'f':
				if ((argv[i][2] != '\0') || (i == argc - 1)) {
					printf ("\nSequence file name not specified or Incorrect usage\n");
					wrong ();
				} else {
					//n_seq=Read_Files(argv[++i],2);
					n_seq = Read_FASTA_File (argv[++i],
							&phi_flag);
					strcpy (fasta_file, argv[i]);
					check[1] = 1;
					if (print_debug) {
						printf ("loading file %s\n",
								fasta_file);
						fflush (stdout);
					}
					break;
				}
			case 'C':
			case 'c':
				codon_counts = 1;
				break;
			case 'P':
			case 'p':
				if ((argv[i][2] != '\0') || (i == argc - 1)) {
					printf ("\nPhi Value file name not specified or Incorrect usage\n");
					wrong ();
				} else {
					Read_Phi_File (argv[++i], n_seq);
					check[17] = 1;
					break;
				}
			case 'U':
			case 'u':{
								 switch (argv[i][2]) {
									 case '1':{
															if ((argv[i][3] != '\0')
																	|| (i == argc - 1)) {
																printf ("\nMutation file name not specified or Incorrect usage\n");
																wrong ();
															} else {
																Read_Mutation_File_1 (argv[++i]);
																strcpy (mutation_file, argv[i]);
																mutation = 1;
																break;
															}
														}

									 case '2':{
															if ((argv[i][3] != '\0')
																	|| (i == argc - 1)) {
																printf ("\nMutation file name not specified or Incorrect usage\n");
																wrong ();
															} else {
																Read_Mutation_File_2 (argv[++i]);
																strcpy (mutation_file, argv[i]);
																mutation = 1;
																break;
															}
														}
									 default:{
														 printf ("Unknown Error... Make sure to specify -U1 or -U2 for incoming sum or relative mutation rate file Exiting...");
														 wrong ();
													 }
								 }
								 break;
							 }
			case 'N':
			case 'n':
							 {
								 switch (argv[i][2]) {
									 case 'e':
									 case 'E':
										 if ((argv[i][3] != '\0') || (i == argc - 1))	//check to make sure it's not followed by another letter 
											 //or it's not the last argument
										 {
											 printf ("\nIncorrect option or no argument given: %s\nExiting\n", argv[i]);
											 wrong ();
										 } else {
											 Ne = atof (argv[++i]);
											 break;
										 }
									 default:
										 printf ("\nIncorrect argument given: %s\nExiting\n", argv[i]);
										 wrong ();
										 break;
								 }
								 break;
							 }
			case 'W':
			case 'w':
							 if ((argv[i][2] != '\0') || (i == argc - 1)) {
								 printf ("\nPrinting file option not specified or Incorrect usage\n");
								 wrong ();
							 } else {
								 check[4] = 1;
								 i++;
								 pout = atoi (argv[i]);
								 if (print_debug)
									 printf ("pout value is %d\n", pout);
								 break;
							 }
			case 'O':
			case 'o':
							 if ((argv[i][2] != '\0') || (i == argc - 1)) {
								 printf ("\nOutput suffix not specified or Incorrect usage\n");
								 wrong ();
							 } else {
								 check[5] = 1;
								 i++;
								 strcpy (out_prefix, argv[i]);
								 break;
							 }
			case 'A':
			case 'a':
							 {
								 switch (argv[i][2]) {
									 case '1':
										 if ((argv[i][3] != '\0')
												 || (i == argc - 1)) {
											 printf ("\nInitiation cost a1 not specified or Incorrect usage\n");
											 wrong ();
										 } else {
											 check[9] = 1;
											 i++;
											 A1 = atof (argv[i]);
											 break;
										 }
									 case '2':
										 if ((argv[i][3] != '\0')
												 || (i == argc - 1)) {
											 printf ("\nElongation cost a2 not specified or Incorrect usage\n");
											 wrong ();
										 } else {
											 check[10] = 1;
											 i++;
											 A2 = atof (argv[i]);
											 break;
										 }
									 case 'T':	//adjust calculations for AT bias. 
										 if ((argv[i][3] != '\0')
												 || (i == argc - 1)) {
											 printf ("\nAT Bias not specified or Incorrect usage. Should be a float\n");
											 wrong ();
										 } else {
											 //check[11]=1;
											 i++;
											 at_bias = atof (argv[i]);
											 if (at_bias >= 1 || at_bias <= 0) {
												 printf ("\nAT Bias %f is out of acceptable range. Must be between 0 and 1.\n", at_bias);
												 wrong ();

											 }
											 break;
										 }
									 default:
										 printf ("\nIncorrect specification of translation costs\n");
										 wrong ();
								 }
								 break;
							 }

			case 'B':
			case 'b':
							 if ((argv[i][2] == 'e' || argv[i][2] == 'E')
									 && (argv[i][3] == 'n' || argv[i][3] == 'N')) {
								 benchmark = 1;
								 break;
							 } else if ((argv[i][2] != '\0')
									 || (i == argc - 1)) {
								 printf ("\nB value not specified or Incorrect usage\n");
								 wrong ();
							 } else {
								 check[11] = 1;
								 i++;
								 B = atof (argv[i]);
								 break;
							 }

			case 'D':
			case 'd':
							 if ((argv[i][2] != '\0') || (i == argc - 1)) {
								 printf ("\nD (print out Delta eta) values not specified or Incorrect usage. Should be 0 or 1\n");
								 wrong ();
							 } else {
								 check[12] = 1;
								 i++;
								 print_delta_eta_vals = atoi (argv[i]);
								 break;
							 }

			case 'M':
			case 'm':
							 //Max time to simulate
							 //If set to a negative value it will be used as a factor times -mu
							 if ((argv[i][2] != '\0') || (i == argc - 1)) {
								 printf ("\nM --max # of time (in generations) to simulate not specified or Incorrect usage.\n");
								 wrong ();
							 } else {
								 check[13] = 1;
								 i++;
								 global_max_time = atof (argv[i]);
								 break;
							 }

			case 'I':	//ignore_aa option. added in version 1.2 
			case 'i':
							 if ((argv[i][2] != '\0') || (i == argc - 1)) {
								 printf ("\nI value not specified or Incorrect usage\n");
								 wrong ();
							 } else {
								 check[14] = 1;
								 i++;
								 ignore_aa = atoi (argv[i]);
								 if (print_debug)
									 printf ("ignoring last %d amino acids\n", ignore_aa);
								 break;
							 }
			case 'R':	//start at random spot if 1
			case 'r':
							 if ((argv[i][2] != '\0') || (i == argc - 1)) {
								 printf ("\nR value not specified or Incorrect usage\n");
								 wrong ();
							 } else {
								 check[15] = 1;
								 i++;
								 random_start = atoi (argv[i]);
								 if (print_debug)
									 printf ("ignoring last %d amino acids\n", ignore_aa);
								 break;
							 }
			case 'V':	//transition to transversion ratio.  
			case 'v':
							 if ((argv[i][2] != '\0') || (i == argc - 1)) {
								 printf ("\nV value not specified or Incorrect usage\n");
								 wrong ();
							 } else {
								 check[15] = 1;
								 i++;
								 gamma_ratio = atof (argv[i]);
								 if (gamma_ratio <= 0) {
									 printf ("V value %f is out of range.  Must be greater than 0\n", gamma_ratio);
									 wrong ();
								 }
								 if (print_debug)
									 printf ("Adjusting transition to transversion ratio to %f\n", gamma_ratio);
								 break;
							 }
			case 'Q':
			case 'q':
							 if ((argv[i][2] != '\0') || (i == argc - 1)) {
								 printf ("\nQ value not specified or Incorrect usage\n");
								 wrong ();
							 } else {
								 check[16] = 1;
								 i++;
								 Q = atof (argv[i]);
								 break;
							 }
			case 'E':
			case 'e':
							 if ((argv[i][2] != '\0') || (i == argc - 1)) {
								 printf ("\nE value not specified or Incorrect usage. Should be 0 or 1 (print eta trace)\n");
								 wrong ();
							 } else {
								 check[18] = 1;
								 i++;
								 print_eta_trace = atoi (argv[i]);
								 //;;printf("value of print_eta_trace: %d\n", print_eta_trace);
								 break;
							 }
		}
	}
}
// Initializing Defaults, Warnings and Errors
for (i = 0; i <= ncheck; i++) {
	if (print_debug) {
		printf ("processing argument %d \n", i);
		fflush (stdout);
	}

	switch (i) {
		case 0:
			if (check[i] == 0) {
				dp = &def_tRNA[0];
				//n_aa=Read_Files(dp,1);
				n_aa = Read_tRNA_File (dp);
				strcpy (trna, def_tRNA);
				break;
			} else if (n_aa < 20) {
				printf ("\n*Warning*\nt_RNA file contains codons for %d amino acids only\n\n", n_aa);
				exit (1);
			}
		case 1:
			if (check[i] == 0) {
				n_seq = Read_FASTA_File (&fasta_file[0],
						&phi_flag);
				if (print_debug) {
					printf ("loading file %s\n",fasta_file);
					fflush (stdout);
				}
			} else if (n_seq < 1) {
				printf ("\nSequence file is empty or not in FASTA format\nCheck README for file formats");
				wrong ();
			}
		case 2:		//carry over from SEMPPR.  Not used here
			if ((check[i] == 0) || (check[i] < 1)) {
				runs = 1000;
				break;
			}
		case 3:
			if (check[i] == 0) {
				MLE = 1;
				break;
			}
		case 4:
			if (check[i] == 0) {
				pout = 1;
				break;
			}
		case 5:
			if (check[i] == 0) {
				strcpy (out_prefix, "output/out");
				break;
			}
		case 6:
			if (check[i] == 0) {
				analytic = 1;
				break;
			}
		case 7:
			if (check[i] == 0) {
				pconf = 0;
				break;
			}
		case 8:
			if (check[i] == 0) {
				gconf = 0;
				break;
			}
		case 9:
			if (check[i] == 0) {
				A1 = 2;
				break;
			}
		case 10:
			if (check[i] == 0) {
				A2 = 4;
				break;
			}
		case 11:
			if (check[i] == 0) {
				B = 0.0025;	//previous to 12/10/08 it was 0.00515
				break;
			}
		case 12:
			if (check[i] == 0) {
				print_delta_eta_vals = 0;
				break;
			}
		case 13:
			if (check[i] == 0) {	//global_max_time = 2E10;
				break;
			}
		case 14:
			if (check[i] == 0) {
				ignore_aa = 0;
				break;
			}
		case 17:
			if (check[i] == 0 && phi_flag == 0) {
				printf ("Phi values not found in fasta file or separate phi file. Exiting...\n");
				wrong ();
			}
	}
}


//set up out_folder based on out_prefix  
j = (int) (strlen (out_prefix)) - 1;
test = 0;
//search backwards for '/'
while (j > 0 && test == 0) {
	if (out_prefix[j] == '/') {
		test = 1;
	} else {
		j--;
	}
}

if (j == 0) {
	//printf("Printing output to local directory\n");
	strcpy (out_folder, "./");
} else {
	for (i = 0; i < j; i++) {
		out_folder[i] = out_prefix[i];
	}
	out_folder[j] = '\0';

	//Check to see if folder exists
	//code from: man 2 stat
	if (stat (out_folder, &sb) == -1) {
		printf ("WARNING: Output folder: \"%s\" does not exist.  Please create.\nExiting\n", out_folder);
		fflush (stdout);
		//perror("stat");
		exit (EXIT_FAILURE);
	}

	if (print_debug) {
		printf ("using output folder: %s\n", out_folder);
		fflush (stdout);
	}
}



//ensure outfolder exists
//from man 3p mkdir.  Make with owner and group r/w 
//mkdir(out_folder, 222); //044 gives d---, 144 gives d-w--w----, 244 gives d-wxrw-r, 344 gives r-x, 224 gives d-wxr-----
//set runs and MLE to 0 if running in analytic mode

if (analytic) {
	runs = 0;
	MLE = 0;
}
//
//  switch(pout)//based on desire output set indicator variables
//    {
//    case -2: //return gmean via stdout
//    case -12: //return id and gmean via stdout
//      calcgmean=1;
//      break;
//    case -11: //return amean via stdout
//    case -1: //return id and amean via stdout
//      calcamean=1;
//      break;
//    case 0: //calc mode-- will be done anyways
//    case -10: //calc mode-- will be done anyways
//      break;
//    default:
//      calcamean = 1;
//      calcgmean = 1;
//      calcvar = 1;
//      break;
//    }      

}


/*======================================================*
 * FUNCTION: Print_Commandline							*
 * 		function to print command line					*
 * 		taken from K&R p. 115							*
 * =====================================================*/

void
Print_Commandline (int argc, char **argv, FILE * outfile) {
	/*while (argc > 0)
		fprintf (outfile, "%s%s", *++argv, (argc > 1) ? " " : "");
	fprintf (outfile, "\n");
	*/

	int i;

	for (i = 0; i < argc; i++) {
		(i == 0) ? fprintf(outfile, "%s", argv[i]) : fprintf(outfile, " %s", argv[i]);
	}
	fprintf(outfile, "\n");
}


/*======================================================*
 * FUNCTION: Print_to_Fasta								*
 * 		Prints evolved Seq to fasta file				*
 * =====================================================*/

void
Print_to_Fasta (struct seq_struct *seq, FILE ** outfile) {
	int i, j;
	int aa_count;
	int max_chr = 70;		//max # of characters/line
	int nt_count;

	aa_count = seq->aa_count;
	fprintf (*outfile, ">%s\tphi=\t%f\n", seq->id, seq->phi_obs);

	nt_count = 0;
	for (i = 0; i <= aa_count; i++) {	//include stop codon
		nt_count += 3;
		if (nt_count < max_chr) {
			fprintf (*outfile, "%s", seq->codon_seq[i]);
		} else {
			nt_count -= 3;
			j = 0;

			while (j < 3) {
				if (nt_count == max_chr) {
					//start new line
					fprintf (*outfile, "\n");
					//reset counter
					nt_count = 0;
				}
				fprintf (*outfile, "%c", seq->codon_seq[i][j]);
				j++;
				nt_count++;
			}
		}
	}
	fprintf (*outfile, "\n");
}




/*======================================================*
 * FUNCTION: Print_Summar_Info							*
 *  	Prints evolution summary info 					*
 * 		(*.evol.summary.tsv)							*
 * =====================================================*/

void
Print_Summary_Info (int id_start, int id_stop, int steps) {
	int i;

	char filename[150];

	FILE *outfile;


	if (id_start < 0) {
		strcpy (filename, out_prefix);
		strcat (filename, ".evol.summary.tsv");
		outfile = fopen (filename, "w+");

		fprintf (outfile, "ORF\t");
		fprintf (outfile, "phi\t");
		fprintf (outfile, "Evolve Steps\t");
		fprintf (outfile, "Evolve Time\t");
		fprintf (outfile, "Eta_initial\t");
		fprintf (outfile, "Eta_final\t");
		fprintf (outfile, "Final \\Delta Eta\t");
		fprintf (outfile, "Mean Mutation \\Delta Eta\t");
		fprintf (outfile, "Var Mutation \\Delta Eta\t");
		fprintf (outfile, "# Fix\n");
		id_start = 0;
		fclose (outfile);
	}

	strcpy (filename, out_prefix);
	strcat (filename, ".evol.summary.tsv");
	outfile = fopen (filename, "a+");

	for (i = id_start; i <= id_stop; i++) {
		fprintf (outfile, "%s\t", Seq[i].id);
		fprintf (outfile, "%.4g\t", Seq[i].phi_obs);
		fprintf (outfile, "%d\t", Seq[i].evol_steps);
		fprintf (outfile, "%.0f\t", Seq[i].evol_time);
		fprintf (outfile, "%.3f\t", Seq[i].eta_initial);
		fprintf (outfile, "%.3f\t", Seq[i].eta_obs);
		fprintf (outfile, "%.3f\t", Seq[i].eta_initial - Seq[i].eta_obs);
		fprintf (outfile, "%.4g\t", Seq[i].delta_eta_mean);
		fprintf (outfile, "%.4g\t", Seq[i].delta_eta_var);
		fprintf (outfile, "\n");
	}

	fclose (outfile);
}


/*======================================================*
 * FUNCTION: Print_Output								*
 * 		Invokes Print_to_Fasta and creates a fasta file *
 * 		and output log									*
 * =====================================================*/

void
Print_Output (struct timeval *time_vec, int argc, char *argv[]) {
	int i, j;
	int time_counter;
	int start_filename;
	int stop_filename;
	FILE *outfile;
	double read_time, comp_time, print_time, total_time;
	char filename[150];
	char tmpstr[60];		//used for testing end of fasta_filename
	//out_prefix is globally defined




	//Calc read time
	time_counter = 1;
	timersub (&time_vec[time_counter], &time_vec[time_counter - 1],
			&time_diff);
	read_time = (time_diff.tv_sec + time_diff.tv_usec / 1000000.);

	//Calc computation time
	timersub (&time_vec[time_counter + 1], &time_vec[time_counter],
			&time_diff);
	comp_time = (time_diff.tv_sec + time_diff.tv_usec / 1000000.);

	time_counter++;

	//print final seq to separate file
	if (pout >= 0) {		//start printout

		//set prefix
		switch (1) {
			case 1:
				strcpy (filename, out_prefix);
				strcat (filename, ".");
				break;

			case 2:
				//put in a separate folder: old approach
				strcpy (filename, out_folder);
				strcat (filename, "/fasta_final/");
				break;

			case 3:		//use a standard name to print output
				strcpy (filename, out_folder);
				i = strlen (filename);
				//check to see if it ends in a /
				//need to offset i by 1 due to C indexing
				if (!(filename[i - 1] == '/'))
					strcat (filename, "/");
				strcat (filename, "ces.output");

		}

		//decide whether to use input file name for creating output file name
		if (0) {
			//process input file name to make output file mimic it in structure
			//get length
			stop_filename = (int) (strlen (fasta_file));

			//get last instance of '/'
			i = stop_filename - 1;

			while (i > 0) {
				if ('/' == fasta_file[i])
					break;
				i--;
			}

			start_filename = ++i;
			for (j = start_filename; j < stop_filename; j++) {
				tmpstr[j - start_filename] = fasta_file[j];
			}
			tmpstr[j - start_filename] = '\0';

			strcat (filename, tmpstr);	//append filename
		}
		//handle fasta at end
		switch (2) {
			case 1:
				//get rid of fasta at end of file name if it exists
				i = (int) (strlen (filename)) - 6;	//move back 6 spots ('fasta/0') 
				for (j = 0; j < 6; j++) {
					tmpstr[j] = filename[i + j];
				}
				tmpstr[j] = '\0';

				//printf("\n\tPrinting to single file. Fasta file ends with %s\n", tmpstr);
				//test to see if fasta_file name ends with 'fasta', if so truncate
				if (strcmp ("fasta", tmpstr)) {
					filename[i] = '\0';	//get rid of '.'
				}
				//append to filename
				strcat (filename, ".fasta");
				break;

			case 2:		//use out prefix and add ces.fasta
				strcat (filename, "ces.fasta");
				break;
			case 3:
				strcat (filename, ".fasta");
				break;
		}

		outfile = fopen (filename, "w+");
		if (outfile == NULL) {
			printf ("Cannot open file %s. Folder must exist before running program. Exiting.\n", filename);
			exit (1);
		}

		fprintf (outfile, "File generated on %s", exec_time);
		fprintf (outfile, "Command Used: \n\t");
		for (i = 0; i < argc; i++) {
			fprintf (outfile, "%s ", argv[i]);
		}
		fprintf (outfile, "\n");

		for (i = 0; i < n_seq; i++) {
			Print_to_Fasta (&Seq[i], &outfile);
		}
		fclose (outfile);
	}



	gettimeofday (&time_vec[time_counter], NULL);
	timersub (&time_vec[time_counter], &time_vec[time_counter - 1],
			&time_diff);
	print_time = (time_diff.tv_sec + time_diff.tv_usec / 1000000.);

	timersub (&time_vec[time_counter], &time_vec[0], &time_diff);
	total_time = (time_diff.tv_sec + time_diff.tv_usec / 1000000.);

	if (print_debug) {
		printf ("Calcuated total_time: %f s\n", total_time);
		fflush (stdout);
	}
	//  timeval_subtract(&print_time, &print_stop,&comp_stop);
	//timeval_subtract(&total_time, &print_stop,&comp_start);
	//  time(&print_stop);
	//print_time=difftime(print_stop,comp_stop);
	//total_time=difftime(print_stop,start);

	//printf("\nCodon Evolution Simulation (CES) Version 0.9\nThe code was executed on %s\nInput file : %s\ntRNA file : %s\nOutput : %s\nNo. of sequences : %d\nNo. of evolutionary steps : %d\nIgnore last AA: %d,\nRuntimes (Read, Computation, Print, Total) : %.4lf\t%.4lf\t%.4lf\t%.4lf seconds\nCommand: %s\n\n\n",exec_time,fasta_file,trna,out_prefix,n_seq, global_max_time,ignore_aa,read_time,comp_time,print_time, total_time, command_line);

	// append log file
	strcpy (filename, out_prefix);
	strcat (filename, ".ces.log");
	outfile = fopen (filename, "a+");
	fprintf (outfile,
			"\nCodon Evolution Simulation (CES) Version 0.9\nThe code was executed on %s\nInput file : %s\ntRNA file : %s\tMutation File : %s\nOutput : %s\nNo. of sequences : %d\nEvolution Run time : %f\nIgnore last AA: %d,\nRuntimes (Read, Computation, Print, Total) : %.4lf\t%.4lf\t%.4lf\t%.4lf seconds\n\n\n",
			exec_time, fasta_file, trna, mutation_file, out_prefix,
			n_seq, global_max_time, ignore_aa, read_time, comp_time,
			print_time, total_time);
	fprintf (outfile, "Command line arguments used to run code: \n\t");
	for (i = 0; i < argc; i++) {
		fprintf (outfile, "%s ", argv[i]);
	}
	fprintf (outfile, "\n");

	fclose (outfile);
	//to get hostname look into hostname.c
	//at  http://www.koders.com/c/fid9A30C5507F7C374C291FEB93CD918B97BEA16C55.aspx?s=md5
	if (print_debug) {
		printf ("Finished printing to log file\n");
		fflush (stdout);
	}

}


/*======================================================*
 * FUNCTION: Evolve_Sequence							*
 * 		This is where the main calculations occur		*
 * =====================================================*/

int
Evolve_Sequence (struct seq_struct *seq, int aa_count) {

	char wt_codon[4];
	char mut_codon[4];

	/* char mut_table_A[3]={'G', 'T', 'C'}; */
	/* char mut_table_T[3]={'C', 'A', 'G'}; */
	/* char mut_table_G[3]={'A', 'T', 'C'}; */
	/* char mut_table_C[3]={'T', 'A', 'G'}; */

	/* char mut_table_a[3]={'g', 't', 'c'}; */
	/* char mut_table_t[3]={'c', 'a', 'g'}; */
	/* char mut_table_g[3]={'a', 't', 'c'}; */
	/* char mut_table_c[3]={'t', 'a', 'g'}; */

	/* char wt_nt; */
	/* char mut_nt; */

	int D_array[MAX_DDIM][2];	//array with codon position and codon index of one step neighbors
	int mut_codon_index;
	int wt_codon_index;
	int mut_codon_pos;	//position numbers for the codon and nt that's mutated


	int D_dim;		//number of dimensions
	int i, j;
	int evol_steps;

	//indices used for Generate_D_Array
	int Dmax;
	int D_index;
	int codon_index;
	//int codon_index_from;
	//int codon_index_to;
	//int aa_index;


	//double Ne = Ne;
	// double mu = MU; //mutation rate per nt
	double phi;		//protein production rate
	double qPhi;
	double two_qNePhi;
	double factor;		//ratio of mut and wt elong_pr, used to rescale sigma_vecs
	double inv_factor;	//inverse of above
	double evol_time;
	double time_step;
	double wait_parameter;
	double max_time;


	//  double sigma_ratio_vec[MAX_AA];
	double b_over_c_vec[MAX_AA];
	double delta_eta_first_term_vec[MAX_AA];
	double delta_eta_second_term_vec[MAX_AA];


	double delta_eta_vec[MAX_DDIM];	//array for storing $\Delta \eta_{i,j}$ values for all
	// one step mutants of the resident allele
	double pi_vec[MAX_DDIM];	//array for storing values of $\pi(i->j)$ for all one step mutants
	double pi_total;		//sum over pi_vec values.
	double mutation_vec[MAX_DDIM];	//array for storing values of mutation rates $\mu_{i->j}$ for all one step mutants
	double pr_vec[MAX_DDIM];	//vector of pi * mu values
	double pr_total;		//sum over pr_vec values.
	double mut_pr;		//RV representing replacement allele
	double curr_pr;		//used to find replacement allele


	//double delta_eta_list[MAX_EVOL_STEPS];
	double delta_eta_mean = 0;	//mean effect of synon sub for the current wt seq-- def differs from before
	double delta_eta_var = 0;	//var effect of synon sub for the current wt seq --definition differs from before
	double delta_eta_sq = 0;	//E(\delta \eta^2), used to calc var


	double evol_delta_eta_mean = 0;	//mean effect of synon sub evolution
	double evol_delta_eta_var = 0;	//var effect of synon sub evolution
	double evol_delta_mean;	//used to update mean and var

	double delta_eta;		//eta_wt-eta_mut

	//double wt_sigma_n;
	double wt_eta;
	//double wt_xi;

	//double mut_sigma_n;
	double mut_eta;
	//double mut_xi;

	//double new_xi,new_sigma_n,new_eta,previous_eta,
	double new_eta,previous_eta,
				 abs_error_delta_eta, rel_error_delta_eta,direct_calc_delta_eta;
	int directly_update_eta = 0;
	//int aa_position,wt_index,mut_index;

	double seq_mu;

	double temp;
	struct timeval curr_time;

	const gsl_rng_type *rng_type;
	gsl_rng *rn;

	char delta_eta_filename[250];
	FILE *delta_eta_outfile;


	char fasta_filename[250];
	FILE *fasta_outfile;

	char eta_trace_filename[250];
	FILE *eta_trace_outfile;

	max_time = seq->max_time;
	//open files for writing output

	//for debugging;
	//int th_id; //thread id
	//th_id = omp_get_thread_num();


	//printf("here with sequence %s\n", seq->id);
	//printf("print_eta_trace is  %d\n", print_eta_trace);
	//printf("print_delta_eta_vals is  %d\n", print_delta_eta_vals);
	//printf("print_evol_fasta is  %d\n", print_evol_fasta);

	fasta_outfile = NULL;

	if (print_evol_fasta) {
		strcpy (fasta_filename, out_folder);
		strcat (fasta_filename, "/fasta_evol");

		//make dir to ensure it exists
		//mkdir(fasta_filename, 664);//S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);

		strcat (fasta_filename, "/");
		strcat (fasta_filename, seq->id);
		strcat (fasta_filename, ".evol.fasta");
		fasta_outfile = fopen (fasta_filename, "w+");
		if(fasta_outfile==NULL){
			printf ("Cannot open file %s. Folder must exist before running program. Exiting.\n", fasta_filename);
			exit (1);
		}
	}


	if (print_delta_eta_vals) {
		strcpy (delta_eta_filename, out_folder);
		strcat (delta_eta_filename, "/delta_eta_lists");

		//make dir to ensure it exists
		//mkdir(delta_eta_filename, 664);//S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);


		strcat (delta_eta_filename, "/");
		strcat (delta_eta_filename, seq->id);
		strcat (delta_eta_filename, ".delta_eta_list.tsv");
		delta_eta_outfile = fopen (delta_eta_filename, "w+");
		if (delta_eta_outfile == NULL) {
			printf ("Cannot open file %s. Folder must exist before running program. Exiting.\n", delta_eta_filename);
			exit (1);
		}
		//print header
		fprintf (delta_eta_outfile,
				"evol_step\tevol_time\tSum_Fix_Pr\teta_wt");
		// Added by Jeremy
		D_dim = Calc_Dimensionality (&(seq->codon_index[0]), aa_count);
		for (i = 0; i < D_dim; i++) {
			fprintf (delta_eta_outfile,
					"\tdelta_eta_%d", (i + 1));
		}
		fprintf (delta_eta_outfile, "\n");
	}

	if(print_eta_trace){
		strcpy (eta_trace_filename, out_folder);
		strcat (eta_trace_filename, "/eta_trace");
		strcat (eta_trace_filename, "/");
		strcat (eta_trace_filename, seq->id);
		strcat (eta_trace_filename, ".tsv");
		//printf("opening eta_trace_filename %s\n", eta_trace_filename);
		eta_trace_outfile = fopen (eta_trace_filename, "w+");
		if (eta_trace_outfile == NULL) {
			printf ("Cannot open file %s. Folder must exist before running program. Exiting.\n", eta_trace_filename);
			exit (1);
		}
		//print header
		fprintf (eta_trace_outfile,
				"ORF\tphi\tevol_step\tdelta_time\teta\tmu_term\n");
	}
	//Set up GSL RNG
	gsl_rng_env_setup ();



	//Create a generator chosen by the 
	//environment variable GSL_RNG_TYPE
	//use default gsl for generating uniform rn 
	//from which other rn functions are derived 

	rng_type = gsl_rng_default;
	rn = gsl_rng_alloc (rng_type);

	//seed rng using clock
	gettimeofday (&curr_time, NULL);
	if (print_debug) {
		printf ("Seeding RNG with clock. Current time is %d %d\n",
				(int) curr_time.tv_sec, (int) curr_time.tv_usec);
		fflush (stdout);
	}

	gsl_rng_set (rn, (curr_time.tv_sec + curr_time.tv_usec));	// (const gsl_rng * r, unsigned long int s)   

	//end RNG set up


	phi = seq->phi_obs;
	qPhi = Q * phi;
	two_qNePhi = qPhi * 2 * Ne;


	//references C++
	//two ways of passing the seq.sigma_vec
	//&seq->sigma_vec
	//&(*seq).sigma_vec
	//However, these seem to pass a pointer to an array as opposed to a pointer
	//to the start of the array, so the syntax that seems to work for my purposes is
	// &(seq->codon_index[0]) gives the exact address I want.

	//seq->value would work for a scalar

	//  sigma_vec = &seq->sigma_vec;

	Convert_Sigma_Vec_to_Sigma_Ratio_Vec (&(seq->sigma_vec[0]),
			&(seq->sigma_ratio_vec[0]),
			aa_count);

	//calculate b_over_c_vec
	Calc_B_over_C_Vec (&(seq->codon_index[0]), b_over_c_vec, aa_count);

	//end Generate_D_Array


	evol_steps = 0;
	time_step = 0;

	for (evol_time = 0; evol_time < max_time; evol_time += time_step) {
		// simulation begins.  Ideally could start further down 
		// and just update below terms in just a few places

		evol_steps++;



		//determining Dimensionality/# one step neighbors of the system.
		//D_dim does change over the simulation!
		//The changes in it's values are due to hte effects of 6 codon aa
		//For example, for Arginine AGA can mutate to AGG or CGA.
		//In contrast, CGA can mutate to AGA, CGC, CGG, or CGU.
		D_dim = Calc_Dimensionality (&(seq->codon_index[0]), aa_count);
		//printf("step = %d; D_dim = %d\n", evol_steps, D_dim);


		//Create D_index of aa position in sequence and synonymous codon index values 
		//Generate_D_Array(seq->codon_index, &D_array, aa_count);
		//  imax = (int)(aa_count); 

		D_index = 0;
		for (i = 0; i < aa_count; i++) {
			codon_index = seq->codon_index[i];
			//aa_index = seq->aa_index[i];
			Dmax = Codon[codon_index].num_one_step_synonym;


			//go through each synonym

			for (j = 0; j < Dmax; j++) {

				mutation_vec[D_index] =
					Codon[codon_index].one_step_relative_mutation_rate[j];

				(D_array[D_index][0]) = (int) (i);//amino acid position
				(D_array[D_index++][1]) =
					(Codon[codon_index].one_step_synonym_index[j]);//synonymous codon for position D_index


			}

		}


		/* if (print_debug && 1 == 1) {
			 printf ("Step: %d: Time of Replacement = %.0f\n\tCalculating Delta_Eta First and Second Terms\n", evol_steps, evol_time + time_step);
			 fflush (stdout);
			 }*/

		Calc_Delta_Eta_NSE_First_and_Second_Term_Vecs (&
				(seq->sigma_ratio_vec
				 [0]),
				b_over_c_vec,
				delta_eta_first_term_vec,
				delta_eta_second_term_vec,
				aa_count);

		if (print_debug && 0 == 1) {
			printf ("\tCalculating Delta Eta Vec\n");
			fflush (stdout);
		}

		Calc_Delta_Eta_NSE_Vec (&(seq->codon_index[0]), b_over_c_vec,
				delta_eta_first_term_vec,
				delta_eta_second_term_vec,
				delta_eta_vec, aa_count, D_dim);

		if (print_debug && 0 == 1) {
			printf ("\tCalculating Pi Vec\n");
			fflush (stdout);
		}

		//Calculate pr of allele replacement given a mutant arises
		Calc_Pi_Vec (delta_eta_vec, pi_vec, &pi_total, phi, qPhi,
				two_qNePhi, aa_count, D_dim);

		//test to see if there are mutation effects
		if ((mu_bias == 0.5) && gamma_ratio == 1. && 0) {
			//flat mutation so pi values = pr values
			//Assign pr_vec address to pi_vec address
			//This should work so long as the dimensions of these
			//vectors are the same
			//      &(pr_vec[0]) = &(pi_vec[0]);//way of copying pointers according to Kernighan and Ritchie p. 98
			//The above doesn't work because pr_vec represents a specific space in memory and it its address can't
			//be reassigned.  It can be copied, however.

			memcpy (pr_vec, pi_vec, D_dim * sizeof (double));

			pr_total = pi_total;

		} else {
			//Calculate true pr of replacement, weighing by mutation rates
			Calc_Replacement_Pr_Vec (mutation_vec, pi_vec, pr_vec,
					&pr_total, D_dim);
		}


		//Calculate expected time to mutation to allele that will replace resident
		//We are assuming that replacement happens quickly.
		wait_parameter = 1 / (Ne * mu * pr_total);	// This is the 'failure' or leaving rate of the resident allele.

		//ideally we would expect the wait parameter to be >>1, but it likely doesn't matter 
		//if(wait_parameter < 100)

		//The time to replacement should follow a geometric dist which we approximate 
		// with an exponential distribution.
		//Note the parameter for GSL's RNG argument is the expected wait time, so it is 1/replacement rate
		time_step = gsl_ran_exponential (rn, wait_parameter);
		if (benchmark||non_random_time_step) time_step = wait_parameter;

		if(print_eta_trace) {
			seq_mu = Calc_Seq_Mu(seq->codon_cts,seq->aa_count);
		}

		/*========================================*
		 *   Print Eta, Mu, and Time Trace Here   *
		 * =======================================*/     

		if(print_eta_trace){
			//print ORF, evolutionary step, time between steps, eta value for current genotype
			fprintf(eta_trace_outfile, "%s\t%g\t%d\t%g\t%g\t%g\n", seq->id, phi, evol_steps, time_step, seq->eta_obs, seq_mu);
		}


		//mut_pr = gsl_rng_uniform(rn) * pr_total;This appears twice- here and line 2534 Whewgley
		if (print_delta_eta_vals) {
			fprintf (delta_eta_outfile, "%d\t%.6e\t%.6e\t%.6e",
					evol_steps, evol_time + time_step, pr_total,
					wt_eta);
			for (i = 0; i < D_dim; i++) {
				fprintf (delta_eta_outfile, "\t%.6e",
						delta_eta_vec[i]);
			}
			fprintf (delta_eta_outfile, "\n");

		}
		//evaluate moments of delta_eta_vec if necessary
		if (print_delta_eta_vals || print_debug) {
			delta_eta_mean = 0.0;
			delta_eta_sq = 0.0;
			for (i = 0; i < D_dim; i++) {
				delta_eta_mean += delta_eta_vec[i];
				delta_eta_sq += pow (delta_eta_vec[i], 2);
			}

			delta_eta_mean /= D_dim;
			delta_eta_sq /= D_dim;

			delta_eta_var = delta_eta_sq - pow (delta_eta_mean, 2);

			seq->delta_eta_mean = delta_eta_mean;
			seq->delta_eta_var = delta_eta_var;
		}else{ 
			//don't calculate this information
			seq->delta_eta_mean = 0;
			seq->delta_eta_var = 0;
		}

		//it seems like for(i=1, i< max_time, i++){ loop could start here if I was to 
		//think about it more.
		// The problem is D_Array set up kinda screwy so that one may need to shift
		// everything up or down depending on the # of neighbors of the mutant
		// Further, lots of stuff needs to be rescaled
		// 
		if (print_debug && 0 == 1) {
			printf ("\tChoosing mutant...");
			fflush (stdout);
		}
		//choose a RV to determine which allele replaces the resident
		mut_pr = gsl_rng_uniform (rn) * pr_total;
		if (benchmark)
			mut_pr = 0.5 * pr_total;

		//Note: pr_total is not scaled by the constants Ne or mu)
		if (print_debug && 0 == 1) {
			printf ("%g from 0 to %g \n", mut_pr, pr_total);
			fflush (stdout);
		}
		//Find out which allele it is by searching forwards. 
		//This is likely more efficient than starting at the end b/c there is weaker selection at 
		// the start of the sequence
		//Some kind of newton search could be even more efficient
		curr_pr = 0.0;
		i = 0;
		while (curr_pr < mut_pr) {
			curr_pr += pr_vec[i++];
		}

		i--;
		D_index = i;


		mut_codon_pos = D_array[D_index][0];
		mut_codon_index = D_array[D_index][1];

		//get wt codon index
		wt_codon_index = seq->codon_index[mut_codon_pos];

		//get codons
		strcpy (wt_codon, Codon[wt_codon_index].codon);
		strcpy (mut_codon, Codon[mut_codon_index].codon);


		delta_eta = delta_eta_vec[i];	//should be \eta_i - \eta_j

		//values to check to make sure delta_eta is working right
		if (1 == 0) {
			wt_eta = Calc_Eta_NSE (seq);
			seq->eta_initial = wt_eta;
		}
		//update everything//

		//update evol_delta_mean/var vals
		evol_delta_mean = delta_eta - delta_eta_mean;
		evol_delta_eta_mean += evol_delta_mean / (evol_time + time_step);
		evol_delta_eta_var +=
			evol_delta_mean * (delta_eta - evol_delta_eta_mean);
		//update values for printing
		seq->evol_delta_eta_mean = evol_delta_eta_mean;
		seq->evol_delta_eta_var = evol_delta_eta_var;

		//update codon 
		strcpy (seq->codon_seq[mut_codon_pos],
				Codon[mut_codon_index].codon);

		//update codon index
		seq->codon_index[mut_codon_pos] = mut_codon_index;

		//update codon counts
		seq->codon_cts[wt_codon_index]--;
		seq->codon_cts[mut_codon_index]++;


		if(directly_update_eta){
			//* WORKAROUND: Calculate eta for new sequence from scratch
			//* this line overwrites the updating of sigma_vec done above which might be the cause of our errors in delta_eta
			//* it should be possible to tell since the delta_eta would work for the first evolutionary step but not the later ones.
			//* TODO: Comment out this line to see if there's a problem with how we update sigma_vec above
			previous_eta= seq->eta_obs;

			new_eta = Calc_Eta_NSE (seq);
			//* this overwrites update of seq->eta_obs using delta_eta
			seq->eta_obs = new_eta;

			if (1 == 1) {
				//* check to make sure delta_eta is working right
				direct_calc_delta_eta = (previous_eta - new_eta);
				abs_error_delta_eta = abs(direct_calc_delta_eta - delta_eta);
				rel_error_delta_eta = abs_error_delta_eta/direct_calc_delta_eta;


				if (rel_error_delta_eta > 1.0E-6) {           
					for(i=0;i<D_dim;i++){
						//aa_position = D_array[i][0];
						//wt_index = seq->codon_index[aa_position];
						//mut_index = D_array[i][1];



					}
					printf("Error: delta eta did not work\n");
				}
			}


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

		//update b_over_c_vec
		b_over_c_vec[mut_codon_pos] =
			B / Codon[mut_codon_index].elong_rate;

		//update eta_obs by subtracting difference: \eta_wt - (\eta_wt - \eta_mut) = \eta_mut
		seq->eta_obs -= delta_eta;


		if (1 == 0) {
			//check to make sure delta_eta is working right
			mut_eta = Calc_Eta_NSE (seq);

			if (abs ((wt_eta - mut_eta) - delta_eta) / delta_eta >
					1.0E-6) {

				printf ("ERROR: wt_eta-mut_eta = %f != %f delta_eta\n", (wt_eta - mut_eta), delta_eta);
			}
		}
		//    eta_obs +=delta_eta;


		if (print_debug || 1 == 0) {
			temp = (evol_time + time_step);
			if (evol_steps == 1)
				printf ("\tevol_time\tmut_codon_pos\twt_codon\tmut_codon\tdelta_eta\tdelta_eta_mean\tdelta_eta_var\n");
			printf ("\t%g\t%d\t%s\t%s\t%g\t%g\t%g\n", temp,
					mut_codon_pos, wt_codon, mut_codon, delta_eta,
					seq->delta_eta_mean, seq->delta_eta_var);


		}
	}				//end sequence evolution



	if (print_evol_fasta) {
		fclose (fasta_outfile);
	}

	if (print_delta_eta_vals) {
		fclose (delta_eta_outfile);
	}

	if (print_eta_trace) {
		fclose (eta_trace_outfile);
	}

	if (pout <= 0) {		//print to individual fasta files for each gene
		//print final seq to separate file
		strcpy (fasta_filename, out_folder);
		strcat (fasta_filename, "/fasta_final/");

		//make dir to ensure it exists
		//mkdir(delta_eta_filename, 664);//S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);

		strcat (delta_eta_filename, "/");

		strcat (fasta_filename, seq->id);
		strcat (fasta_filename, ".final.fasta");

		fasta_outfile = fopen (fasta_filename, "w+");
		if (fasta_outfile == NULL) {
			printf ("Cannot open file %s. Folder must exist before running program. Exiting.\n", fasta_filename);
			exit (1);
		}


		Print_to_Fasta (seq, &fasta_outfile);
		fclose (fasta_outfile);
	}
	//calc mean
	//
	//  delta_eta_mean = 0;
	//  for(i=0; i<num_fix;i++){
	//    delta_eta_mean+=delta_eta_list[i];
	//  }
	//  delta_eta_mean /=num_synon_mut;
	//
	//  //calc var
	//  delta_eta_var = 0;
	//  for(i=0; i<num_fix;i++){
	//    delta_eta_var+=pow(delta_eta_list[i], 2);
	//  }
	//  delta_eta_var /=(num_synon_mut-1);

	//update Seq[]
	seq->evol_steps = evol_steps;
	seq->evol_time = evol_time;
	seq->delta_eta_mean = delta_eta_mean;
	seq->delta_eta_var = delta_eta_var;


	// printf("%s mutations with %d synon and  %d fixations. \nRatios %f and %f\n\n", max_time,num_synon_mut, num_fix, (double)(num_fix)/max_time, (double)(num_fix)/num_synon_mut);
	gsl_rng_free (rn);

	return 1;
	}

	void Codon_Counts(){
		int i;
		int j;
		int total_counts[61] = {0};
		int codon_indx;



		for(i=0;i<n_seq;i++){

			Convert_Codon_Seq_to_Codon_Index(Seq[i].codon_seq,Seq[i].codon_index,Seq[i].aa_count);

			for(j=0;j<Seq[i].aa_count;j++){
				codon_indx = Seq[i].codon_index[j];
				Seq[i].codon_cts[codon_indx]++;
			}
		}

		if(codon_counts) {
			for(i=0;i<n_seq;i++){
				if(Seq[i].phi_obs <= .01){
					for(j=0;j<61;j++){
						total_counts[j] += Seq[i].codon_cts[j];
					}
				}				
			}


			for(i=0;i<61;i++){
				printf("%c\t%s\t%d\n",Codon[i].aa,Codon[i].codon,total_counts[i]);
			}
		}

	}




	/*======================================================*
	 * FUNCTION: main										*
	 * 		Main Function									*
	 * =====================================================*/

	int
		main (int argc, char *argv[]) {
			int i, D_dim;
			int time_counter, aa_count;
			double max_time;
			//double mut_obs;
			time_t start;
			//double time_spent[5]; //time spent in each step
			//struct timeval time_val, time_diff;//, comp_start, comp_stop, print_stop, time_diff;
			struct timeval time_vec[5];	//0 prog_start, 1 read_stop, 2 comp_stop, 3 print_stop
			char tmp_str[100];
			time_t time_raw_format;

			//get run start time
			time (&time_raw_format);
			sprintf (exec_time, "%s", ctime (&time_raw_format));



			// Initializing the random number generator
			srand (time (NULL));	//should this be NULL? mikeg

			time_counter = 0;
			time (&start);
			timeinfo = localtime (&start);

			gettimeofday (&time_vec[time_counter++], NULL);

			//read in any arguments passed on the command line
			Read_Commandline_Args (argc, argv);

			//save command line info to a string
			//based on code from K&R p. 115

			command_line[0] = '\0';
			for (i = 0; i < argc; i++) {
				sprintf (tmp_str, "%s%s", argv[i], (i < argc - 1) ? " " : "");
				strcat (command_line, tmp_str);
			}

			//adjust aa_counts
			if (ignore_aa > 0) {
				for (i = 0; i < n_seq; i++) {
					Seq[i].aa_count -= ignore_aa;
				}
			}

			gettimeofday (&time_vec[time_counter++], NULL);



			Process_AA_Information ();
			//initialize codon structures


			Generate_Codon_Structures ();


			// Main code begns
			Print_Summary_Info (-1, -1, 0);


			if (global_max_time < 0) {
				max_time = -global_max_time / mu * 1.0;

			} else {
				max_time = global_max_time * 1.0;
			}

			//double mutation_vec[MAX_DDIM] = {0};

			Codon_Counts(); //Set initial Codon counts 

			if(!codon_counts) { //Do not run this code if flag is set to only calculate codon counts

				//#pragma omp parallel for private(i, D_dim,aa_count,mutation_vec) shared(Seq,AA,Codon, max_time)
#pragma omp parallel for private(i, D_dim,aa_count) shared(Seq,AA,Codon, max_time)
				//"pragma omp parallel" is the command which will parallelize the next full command (here a for loop)
				// private is defined for individual processes inaccessible by others
				// shared is for variables accessible to all processes


				for (i = 0; i < n_seq; i++) {	//parallelized for loop

					//printf("Starting Sequence %i\n",i);

					aa_count = Seq[i].aa_count;
					//choose outcome
					if (print_debug && 1 == 0) {
						printf ("%d:%s\tProcessing\n", i, Seq[i].id);
						fflush (stdout);
					}


					if (print_debug && 0==1) {
						printf ("%d:%s\tConverting Codon strings to index values\n",
								i, Seq[i].id);
						fflush (stdout);
					}
					Convert_Codon_Seq_to_Codon_Index (Seq[i].codon_seq,
							Seq[i].codon_index, aa_count);

					if (print_debug && 0==1) {
						printf ("Calcuating Dimensionality in main() as %d\n",
								D_dim);
						fflush (stdout);
					}
					//determining Dimensionality/# one step neighbors of the system.
					D_dim = Calc_Dimensionality (Seq[i].codon_index, aa_count);

					if (print_debug && 0==1) {
						printf ("%d:%s\tConverting Codon index values to amino acid index values\n", i, Seq[i].id);
						fflush (stdout);
					}
					Convert_Codon_Index_to_AA_Index (Seq[i].codon_index,
							Seq[i].aa_index, aa_count);

					if (random_start) {
						//replace codn sequence and index  with a randomly generated one
						Generate_Random_Codon_Seq_and_Index (Seq[i].codon_seq,
								Seq[i].codon_index,
								Seq[i].aa_index,
								aa_count);
					}

					//calculate mu_obs for initial sequence
					//To be completed. Need to understand information in AA[i].mu[][] matricies
					//jmax = Seq[i].aa_count;
					//mu_obs = 1;
					//for(j = 1; j<jmax; j++){
					//  codon_index = Seq[i].codon_index[j];
					//  mu_obs*= Codon[codon_index].one_step_relative_mutation_rate
					//}
					//Seq[i].mu_obs = ;


					if (Q > 0.0) {
						////////////////////
						//NSE Model
						////////////////////

						if (print_debug && 0==1) {
							printf ("Using NSE Model with Q = %f, \t\t\tCalculating Sigma Vec\n", Q);
							fflush (stdout);
						}

						Seq[i].eta_initial = Seq[i].eta_obs =
							Calc_Eta_NSE (&(Seq[i]));

						Seq[i].max_time = max_time;

						Evolve_Sequence (&(Seq[i]), aa_count);

					} else {
						////////////////////////////////////////////////
						// ROC Model
						///////////////////////////
						if (print_debug) {
							printf ("Using ROC Model since Q =  %f, \t\t\tCalculating Sigma Vec\n", Q);
							fflush (stdout);
						}

						Seq[i].eta_initial = Seq[i].eta_obs = 
							Calc_Eta_NSE (Seq);

						Seq[i].max_time = max_time;

						Evolve_Sequence (&Seq[i], aa_count);
					}
					Print_Summary_Info (i, i, max_time);

					//      Old_Print_to_Fasta(i);



				}				//end parallel part of code
			}

			if (print_debug) {
				printf ("Finished calculations... printing output..\n");
				fflush (stdout);
			}


			gettimeofday (&time_vec[time_counter++], NULL);

			if(!codon_counts){
				Print_Output (time_vec, argc, argv);
			}


			return 0;


		}
