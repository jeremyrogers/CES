#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "CES.h"

int Read_FASTA_File (char *filename, int *phi_flag);
extern int Read_tRNA_File(char *filename);// located in CES.h??
void Read_Mutation_File_1 (char *filename);

struct amino_acid AA[22];
struct codon_struct Codon[64];
struct seq_struct Sequence;

int main(){

char fasta_file[100], trna_file[100], mut_file[100];
int * phi_flag;
//define filenames
sprintf(fasta_file,"test.fasta");
sprintf(trna_file,"test.tRNA");
sprintf(mut_file,"test.mut");

//read in files
Read_FASTA_File("test.fasta",phi_flag);
//Read_tRNA_File(trna_file);





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
      char *nt_seq;//[10000];
      nt_seq = (char *) malloc(sizeof(char)*10000);
      int end_of_seq;
      int nt_count, aa_count;
      FILE *fh;

      fh = fopen (filename, "r");
      if (!fh) {
	    printf ("\nFASTA/Sequence File Doesn't Exist\n");
	    //wrong ();
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
	    while ((curr_char != ' ') && (curr_char != '\t')
		   && (curr_char != '\n')
		   && (k < 26)) {
		  curr_char = curr_line[k + 1];	//offset by 1 b/c of >
		  Sequence.id[k] = curr_char;
		  k++;
	    }
	    if (k == 25) {
		  printf ("Seq ID greater than 25 characters. Exiting\n");
		  exit (1);
	    }
	    //replace new line with null character
	    Sequence.id[--k] = '\0';

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
		  Sequence.phi_obs = 0;//strtod (char_ptr, NULL);	//first argument is ptr to str location of double
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
				printf ("Seq[%d].id = %s has %d nts.  MAX_NTS is %d. Exiting", i, Sequence.id, nt_count, MAX_NTS);
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
		  printf ("Seq[%d].id = %s had %d nts, which is not a multiple of 3. Exiting", i, Sequence.id, nt_count);
		  exit (1);
	    } else {		//everything seems okay
		  aa_count = (nt_count / 3) - 1;	//subtract 1 b/c of stop codon
		  Sequence.aa_count = aa_count;

		  k = 0;

		  for (j = 0; j <= aa_count; j++) {	//include stop codon
			Sequence.codon_seq[j][0] = nt_seq[k++];
			Sequence.codon_seq[j][1] = nt_seq[k++];
			Sequence.codon_seq[j][2] = nt_seq[k++];
			Sequence.codon_seq[j][3] = '\0';
		  }

	    }

	    i++;		//increment seq id

      }
      free(nt_seq);
      return i;
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
      int i, j, k,l;
      FILE *file_handle;
      char curr_char;
      char curr_codon[4];
      double mu_rate;
      int num_codons;
      int max_aa = 22;		//20 AA but serine may get split and we may have stop codons
      int codon_index;
      int test;

      j = 0;

      file_handle = fopen (filename, "r");
      if (!file_handle) {
	    printf ("\nMutation File: %s Doesn't Exist\n", filename);
	    //wrong ();
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
			
		  fscanf (file_handle, "%lf", &mu_rate);

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
			printf ("\nCurrent codon %s did not match any codons for aa %s. Exiting...\n", curr_codon, AA[i].aa);
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
