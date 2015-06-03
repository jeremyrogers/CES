#ifndef IOFUNS_CES_Header_
#define IOFUNS_CES_Header_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
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
 
 int Read_tRNA_File (char *filename);
 int Process_AA_Information ();
 void Generate_Codon_Structures (); 
 int wrong();
 
#endif /* IOFUNS_CES_Header_ */
