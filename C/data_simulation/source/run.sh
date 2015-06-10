#! /bin/bash

################################################################################
#                                                                              #  
# MAKE SURE THIS SCRIPT IS IN THE SAME WORKING DIRECTORY AS THE CES EXECUTABLE #
#                                                                              #
################################################################################

dir=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )


###################################################################
#                                                                 #
#  File locations: give empty string to set to default or ignore  #
#                                                                 #
###################################################################

datadir="../../../data"

# fasta file
genome="$datadir/fasta/jeremyYeast.fasta"

# tRNA abundance/codon translation file
trna="$datadir/elong_rates/jeremyYeast.splitSer.tsv"

# phi value file, unnecessary if phis in fasta
phi="$datadir/expression_rates/jyeast.phi.tsv"

# Output codon counts? (true/false)
ccounts=false

# incoming sum mutation rate file
summut=""

# individual mutation rate file
indmut=""

# ribosome initiation cost in ATPs
# A_1 in SEMPPR (DEFAULT = 4)
A_1=""

# ribosome elongation cost in ATPs
# A_2 in SEMPPR (DEFAULT = 4)
A_2=""

# genome AT bias (0.0-1.0) (DEFAULT = 0.5)
AT=""

# B parameter value (DEFAULT = 0.0025)
B=""

# Print delta_eta files for each evolutionary step 
# Print outs can get quite large (true/false)
deltaeta=false

# Print out eta and mu trace files (true/false)
etamu=false

# CES will relax selection on the last AAs of each sequence
# Specify how many
relax=""

# Simulation time, scaled by mutation rate mu (DEFAULT = -20)
simtime=""

# Effective population size (DEFAULT = 1.36E7)
popsize=""

# output files (DEFAULT = ./output/out)
output=""

# Ratio of transition to transversion mutation rates (DEFAULT = 1.0)
ratio=""

# Replace randomly generated numberes with predetermined ones (true/false)
ben=false

syscall="$dir/CES"

if [ ! -z "$genome" ] ; then
	syscall="$syscall -F $genome"
fi

if [ ! -z "$trna" ] ; then
	syscall="$syscall -T $trna"
fi

if [ ! -z "$phi" ] ; then
	syscall="$syscall -P $phi"
fi

if [ "$ccounts" = true ] ; then
	syscall="$syscall -C"
fi

if [ ! -z "$summut" ] ; then
	syscall="$syscall -U1 $summut"
fi

if [ ! -z "$indmut" ] ; then
	syscall="$syscall -U2 $indmut"
fi

if [ ! -z "$A_1" ] ; then
	syscall="$syscall -A1 $A_1"
fi

if [ ! -z "$A_2" ] ; then
	syscall="$syscall -A2 $A_2"
fi

if [ ! -z "$AT" ] ; then
	syscall="$syscall -AT $AT"
fi

if [ ! -z "$B" ] ; then
	syscall="$syscall -B $B"
fi

if [ "$deltaeta" = true ] ; then
	syscall="$syscall -D 1"
fi

if [ "$etamu" = true ] ; then
	syscall="$syscall -E 1"
fi

if [ ! -z "$relax" ] ; then
	syscall="$syscall -I $relax"
fi

if [ ! -z "$simtime" ] ; then
	syscall="$syscall -M $simtime"
fi

if [ ! -z "$Ne" ] ; then
	syscall="$syscall -Ne $Ne"
fi

if [ ! -z "$output" ] ; then
	syscall="$syscall -O $output"
fi

if [ ! -z "$ratio" ] ; then
	syscall="$syscall -V $ratio"
fi

if [ "$ben" = true ] ; then
	syscall="$syscall -ben"
fi

echo "$syscall"

$syscall
