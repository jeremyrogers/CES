#include "iofuns_CES.h"

/*=======================================================================*
 * FUNCTION:Read_tRNA_File                                               *
 *      Purpose: This function reads a tRNA file and stores information  *
 *               in amino_acid and codon_struct structures.              *
 *               See readme for details concerning file formats          *
 *      Inputs:  char *filename - string containing tRNA filename        *
 *      Outputs: int aa_processed - number of amino acids in tRNA file   *
 *      Usage:   sprinf(filename,"path/to/tRNA_file");                   *
 *               num_aa = Read_tRNA_File(filename);                      *
 * ======================================================================*/
 
int Read_tRNA_File (char *filename) 
{
      int i, j;
      FILE *file_handle;
      char curr_char;
      char aa_list[] = { 'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'Z', 'S', 'T', 'V', 'W', 'Y', 'X' };    //'X' represents the 'stop' amino acid, 'Z' represents the smaller set of serine codons: AGT and AGC.  The other ser codons are TCT, TCC, TCG, and TCA.
      //Note program expects 'X' to be at end of list so 'Z' is put in internally.
      int num_codons;
      int aa_processed;
      int max_aa = 22;      //20 AA but serine may get split and we may have stop codons
			int fscanf_check;



      j = 0;
      //initialize AA structures
      for (i = 0; i < max_aa; i++) {
        AA[i].aa = aa_list[i];
        AA[i].num_codons = 0;
      }

      file_handle = fopen (filename, "r");
      if (!file_handle) {
        if(R_or_C == 1){
          printf ("\ntRNA File: %s Doesn't Exist\ntRNA", filename);
          wrong ();
        }else{
          Rprintf("\ntRNA File: %s Doesn't Exist\ntRNA", filename);
          exit(EXIT_FAILURE);
        }
        
      }
      //get aa index 
      curr_char = fgetc(file_handle);

      while (curr_char != EOF) {
        i = 0;      //set AA index to 0
        //flip through AA until you get a match
        while (curr_char != AA[i].aa && i < max_aa) {
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
            if(R_or_C == 1){
              printf ("\nAA index '%c' in file did not match any known AA. Exiting...\n",curr_char);
            }else{
              Rprintf ("\nAA index '%c' in file did not match any known AA. Exiting...\n",curr_char);
            }
            exit (1);
          }
        } else {
          //get current codon count
          num_codons = AA[i].num_codons;

          if (num_codons >= 6)  //check num_codons value 
          {
            printf ("\nCodon count for AA is greater than maximum possible value of 6.  Exiting...\n");
            exit (1);
          }


          curr_char = fgetc (file_handle);  //get next character. Should be a \t

          //read in codon
          for (j = 0; j < 3; j++)   //load codon sequence and put in codon_index codon index
          {
            AA[i].codon[num_codons][j] = fgetc (file_handle);
          }
          AA[i].codon[num_codons][j] = '\0';


          curr_char = fgetc (file_handle);  //get next char. Should be a \t

          //read in elongation rate
          fscanf_check = fscanf (file_handle, "%lf", &AA[i].elong_rate[num_codons]);

					if (fscanf_check != 1) {
						fprintf(stderr, "Reading in elongation file failed\n");
						exit(1);
					}


          //elong_rate = AA[i].elong_rate[num_codons];
          //increment codon count for the aa
          AA[i].num_codons++;
          num_codons++;


          //read until end of line or EOF
          while ((curr_char != '\n') && (curr_char != EOF))
            curr_char = fgetc (file_handle);

          //get aa index for next codon if not at EOF
          if (curr_char != EOF)
            curr_char = fgetc (file_handle);
        }
      }
      fclose (file_handle);


      //check to make sure the correct # of AA and codons are defined
      aa_processed = 0;
      for (i = 0; i < max_aa; i++) {
        if (AA[i].num_codons > 0) {
          aa_processed++;
        }
      }
      return aa_processed;  //return # of 

}

/*==============================================================*
 * FUNCTION: Process_AA_Information                             *
 *      Purpose: Find max/min elongation rates, calc elong_pr,  *
 *               fail_pr                                        *
 *      Inputs: None                                            *
 *      Outputs: return(1) if the fcn successfully exits        *
 *      Usage: Process_AA_Information();                        *
 * =============================================================*/
 
int Process_AA_Information () {
      char nt;
      char codon[4];
      int i, j, k;
      int num_codons;
      int max_aa = 22;
      double mu_in = 0;
      double sum_neutral_obs_pr;
      double elong_rate, min_elong_rate, max_elong_rate;

      mu_bias = at_bias / (1 - at_bias);    //calculate biased mutation rate based on AT bias
      //Added 3/2/09
      //check to make sure the correct # of AA and codons are defined
      j = 0;
      k = 0;

      if (print_debug) {
         if (R_or_C == 1) {
            printf ("Processing AA information\n");
            fflush (stdout);
         } else {
            Rprintf ("Processing AA information\n");
         }
      }
      //check for if stop codons defined. If not, define them.
      if (AA[max_aa - 1].num_codons == 0) {
        AA[max_aa - 1].num_codons = 3;
        strcpy (AA[max_aa - 1].codon[0], "TGA");
        strcpy (AA[max_aa - 1].codon[1], "TAG");
        strcpy (AA[max_aa - 1].codon[2], "TAA");
      }

      for (i = 0; i < max_aa; i++) {
        k += AA[i].num_codons;
      }

      if (k != 64) {
          if (R_or_C ==1) {
            printf ("Only %d codons defined. Expecting 64. Exiting", k);
            exit (1);
          } else {
            Rprintf ("Only %d codons defined. Expecting 64. Exiting", k);
            exit (1);
          }
      }

      for (i = 0; i < max_aa - 1; i++) {
        num_codons = AA[i].num_codons;

        sum_neutral_obs_pr = 0; //added line 11/21/08

        for (j = 0; j < num_codons; j++) {
          elong_rate = AA[i].elong_rate[j];
          AA[i].elong_pr[j] = elong_rate / (elong_rate + B);
          if (print_debug && 0==1) {
              if (R_or_C){
                printf ("%d:%d %s: %f\n", i, j, AA[i].codon[j], AA[i].elong_pr[j]);
                fflush (stdout);
              } else {
                Rprintf ("%d:%d %s: %f\n", i, j, AA[i].codon[j], AA[i].elong_pr[j]);
              }
          }
          AA[i].fail_pr[j] = B / (elong_rate + B);
          //first part of calculating the pr of obs codon under neutral model but with AT bias
          //set to 1 and then scale based on at_bias
          AA[i].neutral_obs_pr[j] = 1;
	    }
	  }
      
      /*=============================================*
       * Function Call: Generate_Codon_Structures()  *
       * We need to generate these structures before *
       * calculating neutral observation pr b/c      *
       * they hold relevant mutation information     *
       *=============================================*/

      Generate_Codon_Structures();
      
      for (i = 0; i < max_aa - 1; i++) {
        num_codons = AA[i].num_codons;

        sum_neutral_obs_pr = 0; //added line 11/21/08

        for (j = 0; j < num_codons; j++) {

          /* More complex Mutation model if mutation file */
          if (mutation) {
            
            if(AA[i].num_codons == 1){
                AA[i].neutral_obs_pr[j] = 1;
            }else{
                mu_in = 0;
                for (k = 0; k < num_codons; k++) {
                      mu_in += AA[i].mu[k][j];
                }
                AA[i].neutral_obs_pr[j] = mu_in;
            }
          } else {
            strcpy (codon, AA[i].codon[j]);
            for (k = 0; k < 3; k++) {
                  nt = codon[k];
                  if (nt == 'A' || nt == 'T') {
                    AA[i].neutral_obs_pr[j] *= mu_bias; // old, wrong code: at_bias;
                  } //else{
                  //AA[i].neutral_obs_pr[j] *= (1-at_bias);
                  //}
            }
          }

          //keep track of total for an AA.  This will be used to scale previous calculations
          sum_neutral_obs_pr += AA[i].neutral_obs_pr[j];
          AA[i].neutral_obs_pr_cum[j] = sum_neutral_obs_pr;
        }

        //scale each codon by sum_neutral_obs_pr to ensure they sum to 1
        for (j = 0; j < num_codons; j++) {
          AA[i].neutral_obs_pr[j] /= sum_neutral_obs_pr;
          AA[i].neutral_obs_pr_cum[j] /= sum_neutral_obs_pr;
        }
      }


      if (print_debug) {
          if (R_or_C) {
            printf ("\tFinding min/max elongation rates and pr\n");
            fflush (stdout);
          } else {
            Rprintf ("\tFinding min/max elongation rates and pr\n");
          }
      }
      
      //Assign min/max_elong_rate/pr
      for (i = 0; i < max_aa - 1; i++) {
        num_codons = AA[i].num_codons;

        min_elong_rate = AA[i].elong_rate[0];
        max_elong_rate = min_elong_rate;


        //start at 1 b/c already processed first codon
        for (j = 1; j < num_codons; j++) {

          if (min_elong_rate > AA[i].elong_rate[j]) {
            min_elong_rate = AA[i].elong_rate[j];
          }

          if (max_elong_rate < AA[i].elong_rate[j]) {
            max_elong_rate = AA[i].elong_rate[j];
          }

        }

        AA[i].min_elong_rate = min_elong_rate;
        AA[i].max_elong_rate = max_elong_rate;

        AA[i].min_elong_pr = min_elong_rate / (min_elong_rate + B);
        AA[i].max_elong_pr = max_elong_rate / (max_elong_rate + B);
      }

// Calculating the first and second moments for each amino acid
      if (print_debug) {
          if (R_or_C) {
            printf ("calculating moments (expected and var in elongation pr)\n");
            fflush (stdout);
          } else {
            Rprintf ("calculating moments (expected and var in elongation pr)\n");
          }
      }
      for (i = 0; i < max_aa - 1; i++) {
        if (print_debug) {
            if (R_or_C){
              printf ("calculating values for aa %d\n", i);
              fflush (stdout);
            } else {
              Rprintf ("calculating values for aa %d\n", i);
            }
        }
        AA[i].e_invc = 0;
        AA[i].e_cpbinvc = 0;
        AA[i].e_invc2 = 0;
        AA[i].e_cpb2invc2 = 0;
        AA[i].e_cpbinvc2 = 0;
        for (j = 0; j < AA[i].num_codons; j++) {
          AA[i].e_invc +=
            1 / (AA[i].elong_rate[j]) * AA[i].neutral_obs_pr[j];
          AA[i].e_cpbinvc +=
            (B +
             AA[i].elong_rate[j]) /
            (AA[i].elong_rate[j]) * AA[i].neutral_obs_pr[j];
          AA[i].e_invc2 +=
            1 / (pow (AA[i].elong_rate[j], 2)) *
            AA[i].neutral_obs_pr[j];
          AA[i].e_cpb2invc2 +=
            (pow (B + AA[i].elong_rate[j], 2)) /
            (pow (AA[i].elong_rate[j], 2)) *
            AA[i].neutral_obs_pr[j];
          AA[i].e_cpbinvc2 +=
            (B +
             AA[i].elong_rate[j]) /
            (pow (AA[i].elong_rate[j], 2)) *
            AA[i].neutral_obs_pr[j];
        }
      }

      if (print_debug) {
          if (R_or_C){
            printf ("exiting moments\n");
            fflush (stdout);
          } else {
            Rprintf ("exiting moments\n");
          }
      }

      return (1);       //exited successfully
}


/*==============================================================*
 * FUNCTION: Generate_Codon_Structures                          *
 *      Purpose: Store aa_index,elong_rate/pr, and info about   *
 *               relative mutation in codon structure--         *
 *      Inputs: None                                            *
 *      Outputs: return(1) if fcn exits correctly               *
 * =============================================================*/

void Generate_Codon_Structures()
{
    int num_syn, num_one_step_syn;
    int num_two_step_synonym;
    double tsn_penalty = 100; //two step neighbor mutation rate will be smaller than osn by this factor
    int i, j, k, l, m;
    int num_matches;
    int codon_index, synon_index, max_aa = 22;
    int c_index;
    int codon_match[3]; //use to determine which position and nts involved in differences b/w codons.
    double relative_mut_rate, mu_rate;
    int to_codon,to_codon_matcher,from_codon;
    int osn_flag;
    int first_codon;
    int from_matrix_index, to_matrix_index;

    char codon[4];
    char synonym[4];
    char orig_nt, mut_nt;

    codon_index = 0;

    for (i = 0; i < max_aa; i++) {
        for (j = 0; j < AA[i].num_codons; j++) {
            AA[i].codon_index[j] = codon_index;
            strcpy(Codon[codon_index].codon, AA[i].codon[j]);
            Codon[codon_index].aa = AA[i].aa;
            Codon[codon_index].aa_index = i;
            Codon[codon_index].incoming_mu = AA[i].incoming_mu[j];
            Codon[codon_index].elong_rate = AA[i].elong_rate[j];
            Codon[codon_index].elong_pr = AA[i].elong_pr[j];
            Codon[codon_index++].fail_pr = AA[i].fail_pr[j];
        }
    }

    //Define synonyms
    for (i = 0; i < 64; i++) {
        num_syn = 0;
        for (j = 0; j < 64; j++) {
            if (j != i) {
                if (Codon[i].aa_index == Codon[j].aa_index) {
                    Codon[i].synonym_index[num_syn++] = j;
                }
            }
        }
        Codon[i].num_synonym = num_syn;
    }

    //Define one step synonyms
    for (i = 0; i < 64; i++) {
        num_syn = Codon[i].num_synonym;
        strcpy(codon, Codon[i].codon);
        num_one_step_syn = 0;

        //go through each syn
        for (j = 0; j < num_syn; j++) {
            num_matches = 0;

            synon_index = Codon[i].synonym_index[j];
            strcpy(synonym, Codon[synon_index].codon);
            //go through each nt in codon
            for (k = 0; k < 3; k++) {
                if (codon[k] == synonym[k]) {
                    num_matches++;
                    codon_match[k] = 1;
                } else {
                    codon_match[k] = 0;
                }
            }
            if (num_matches == 2) {
                Codon[i].one_step_synonym_index
                    [num_one_step_syn] = synon_index;

                //find mismatch
                k = 0;
                while (codon_match[k] == 1) {
                    k++;
                }

                //get nts
                orig_nt = codon[k];
                mut_nt = synonym[k];

                //mikeg: following code does two things: 
                // 1) Calculate mtuation rates under simple 
                //     AT bias and Transition vs. Transversion bias
                // 2) Assign relative mutation values for one step neighbors
                //Function 1) should be put into a separate function that's
                // used when no mutation file is read in and calculations should
                // should be stored in AA.mu[][] structures as is done when reading
                // mutation information in from a file
                //Function 2) should be done here and it should use the info 
                // from AA.mu[][]
                //
                //As it stands, it appears that the information that Preston
                // assigned to the Preston's code reads in the mu files
                // writes this info to the AA.mu[][], but does not use this
                // information here
                //

                //set initial relative mutation rate

                relative_mut_rate = 1;

                switch (mut_nt) {
                case 'A':
                case 'a':
                case 'T':
                case 't':
                    relative_mut_rate *= mu_bias;
                    break;

                default:
                    //Do nothing
                    break;
                    //relative_mut_rate*=(1-at_bias);
                }   //end case

                //determine if rate is scaled by transition/transversion ratio: gamma_ratio
                if ((orig_nt == 'G' && mut_nt == 'A')
                    || (orig_nt == 'A' && mut_nt == 'G')
                    || (orig_nt == 'C' && mut_nt == 'T')
                    || (orig_nt == 'T' && mut_nt == 'C')) {
                    relative_mut_rate *= gamma_ratio;
                }

                Codon[i].one_step_relative_mutation_rate
                    [num_one_step_syn] = relative_mut_rate;

                num_one_step_syn++; //increase count

            }   //end match for one step synonym

        }
        Codon[i].num_one_step_synonym = num_one_step_syn;
    }

    if (mutation) {
        
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
        
        //Initialize mu matrix in AA struct
        for(i=0;i<max_aa;i++){
            for(j=0;j<AA[i].num_codons;j++){
                for(l=0;l<AA[i].num_codons;l++){
                    AA[i].mu[j][l]=0;
                }
            }
        }
        for (i = 0; i < max_aa; i++) {
            
            first_codon = AA[i].codon_index[0];
            
            for (j = 0; j < AA[i].num_codons; j++) {
				
				c_index = AA[i].codon_index[j];
				num_syn = Codon[c_index].num_synonym;
				num_one_step_syn = Codon[c_index].num_one_step_synonym;
				to_codon = c_index;
			
				mu_rate = 0;
	
				if(AA[i].num_codons > 1){
					if(num_syn==num_one_step_syn){
						//This is the mu rate FROM one-step synonyms TO Codon[c_index]
						mu_rate = Codon[c_index].incoming_mu/num_one_step_syn;
						to_codon = c_index;
					}else{
						num_two_step_synonym = num_syn - num_one_step_syn;
						
						//This is the mu rate FROM one-step synonyms TO Codon[c_index]
						//mu rate FROM two-step synonyms TO Codon[c_index] = mu_rate/tsn_penalty
						mu_rate = Codon[c_index].incoming_mu/(num_one_step_syn + num_two_step_synonym/tsn_penalty);
						to_codon = c_index;
					}
				}
				
				from_codon = 0; // to make the compiler happy

				for(k = 0; k < Codon[c_index].num_synonym;k++){
					l=0;
					
					from_codon = Codon[c_index].synonym_index[k];
					
					//* Find out if codon is a one step neighbor *//
					osn_flag = 0;
                    for(m=0;m<Codon[c_index].num_one_step_synonym;m++){
						if(Codon[c_index].one_step_synonym_index[m] == from_codon) osn_flag = 1;
					}
									
					to_codon_matcher = 0; //* This is what we will use to search for the "to" codon in the "from" codon structure
                   
                    do{ to_codon_matcher = Codon[from_codon].synonym_index[l++];}while(to_codon_matcher != to_codon);
                    l--; //counteract l++ in while loop
                                        
                    //* from_codon = Codon[from_codon]
                    //* to_codon = Codon[to_codon]
                    //* from_matrix_index = from_codon minus first_to_codon
                    //* to_matrix_index = codon_index minus first_codon_index
                    
                    from_matrix_index = from_codon - first_codon;
                    to_matrix_index = to_codon - first_codon;
                    
                    if(osn_flag){                    
						Codon[from_codon].syn_relative_mutation_rate[l] = mu_rate; 
						AA[i].mu[from_matrix_index][to_matrix_index] = mu_rate;
					}else{
						Codon[from_codon].syn_relative_mutation_rate[l] = mu_rate/tsn_penalty;
						AA[i].mu[from_matrix_index][to_matrix_index] = mu_rate/tsn_penalty;
					}
					
                    
					
				}
                
                
                

            }
        }
        
        /*for(i=0;i<max_aa;i++){
            num_codons = AA[i].num_codons;
            for (j = 0; j < num_codons; j++) {
                codon_index = AA[i].codon_index[j];
                for(k=0;k<Codon[codon_index].num_one_step_synonym;k++){
                    Rprintf("%f\t",Codon[codon_index].one_step_relative_mutation_rate[k]);
                }
                Rprintf("\n");
            }
        }*/
        
        //Import from AA.mu matrix to Codon.one_step_relative_mutation vector 
        /*for (i = 0; i < max_aa; i++) {
            num_codons = AA[i].num_codons;
            for (j = 0; j < num_codons; j++) {
                counter = 0;
                codon_index = AA[i].codon_index[j];
                for(l=0;l<Codon[codon_index].num_one_step_synonym;l ++) {
                    for (k = 0; k < num_codons; k++) {
                        if (Codon[j].one_step_synonym_index[l] == AA[i].codon_index[k] ) {
                            Codon[codon_index].
                                one_step_relative_mutation_rate
                                [counter++] =
                                AA[i].mu[j][k];
                        }
                    }
                }
            }
        }*/
        //Check to see if it worked
/*      for(i=0;i<max_aa;i++){
            num_codons = AA[i].num_codons;
            for(l=0;l<num_codons;l++){
                Rprintf("\n");
                for(k=0;k<num_codons;k++){
                    Rprintf("%f\t",AA[codon_index].mu[l][k]);
                }
            }
            Rprintf("\n");
        }*/
        
    
    }
    //Define factors of (c_i/(c_i+b))/(c_j/(c_j+b)) for each one step synonym
    for (i = 0; i < 64; i++) {
        num_syn = Codon[i].num_synonym;
        //go through each one step syn
        for (j = 0; j < num_syn; j++) {
            synon_index = Codon[i].synonym_index[j];
            Codon[i].elong_pr_ratio[j] =
                Codon[i].elong_pr / Codon[synon_index].elong_pr;
        }
    }

    //Define differences of (b/c_i - b/c_j) for each one step synonym
    for (i = 0; i < 64; i++) {
        num_syn = Codon[i].num_synonym;
        //go through each one step syn
        for (j = 0; j < num_syn; j++) {
            synon_index = Codon[i].synonym_index[j];
            Codon[i].delta_b_over_c[j] =
                B * (1. / Codon[i].elong_rate -
                 1. / Codon[synon_index].elong_rate);
        }
    }

}

/*========================================================*
 * FUNCTION: wrong()                                      *
 *      Purpose: This function prints an error message    *
 *               and exits the program. Used in           *
 *               Read_Commandline_Args and Read_tRNA_File *
 *      Inputs: NONE                                      *
 *      Outputs: NONE                                     *
 *      Usage: printf("Error...");                        *
 *             wrong();                                   *
 *========================================================*/
 
int wrong () {
    if(R_or_C == 1){
      printf ("\nError in command line:\nCES Exiting...\n");
    }else{
      Rprintf("\nError in command line:\nCES Exiting...\n");
    }
    exit(1);
}
