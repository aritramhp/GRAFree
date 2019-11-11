# sample script file
# author - Aritra Mahapatra
# date: 22.03.2019

# Modification-3:
#			1) Segment the point set of drift into m number of segments, P[0]...P[m-1].
#			2) Compute the features vector for each segment 
#			3) Compute the distance using the same metric for P[i] and Q[i] 
#			4) Final distance will be the sum of the distances of the segments.
#			5) Make a list of fragments, thus this program can run for multiple fragments 

MODE_1=1				# MODE_1: Compute the composition of genome using Genome_Composition.py

MODE_2=1				# MODE_2: 1 = run GRAFree.py
					
MODE_3=1				# MODE_3: 1 = Compute the bootstrap support scores using SumTree.py


#------------------------------------------------------------------------------------------------
# Declaration of BLOCK and fragments
#------------------------------------------------------------------------------------------------
BLOCK_SIZE=550      # Declaration of size of BLOCK
NO_OF_FRAGMENTS=15    # Declaration of number of fragments

#-------------------------------------------------------------------------------------------------
# Declaration of the BASE Locations of different directories
#-------------------------------------------------------------------------------------------------
DATASET_BASE_DIR='Sample_Execution'
OUTPUT_BASE_DIR='Sample_Execution/Output'

#--------------------------------------------------------------------------------------------------
# Declaration of number of replica of bootstrap
#--------------------------------------------------------------------------------------------------
num_of_replica=10

#--------------------------------------------------------------------------------------------------
# Declaration of save the plot in this format
#--------------------------------------------------------------------------------------------------
SAVE_DATAFORMAT='pdf'

#--------------------------------------------------------------------------------------------------
# Declaration of DATASET
#--------------------------------------------------------------------------------------------------
DATASET='sample_dataset'
DATASET_DIR=$DATASET_BASE_DIR'/'$DATASET
OUTPUT_DIR=$OUTPUT_BASE_DIR'/'$DATASET
INPUT_FILENAME='complete_sequence.fasta'

#--------------------------------------------------------------------------------------------------
# Declaration of ASTRAL PACKAGE.
#--------------------------------------------------------------------------------------------------
ASTRAL_PKG='ASTRAL-master-4.10.12/astral.4.10.12.jar'

#-----------------------------------------------------------------------------------
# This script compute the composition of genome
#-----------------------------------------------------------------------------------
if [ $MODE_1 == 1 ]
then
	GENOME_COMP_OUTDIR=$OUTPUT_DIR'/Genome_Composition'
	execstr='./Genome_Composition.py -D '$DATASET_DIR' -I '$INPUT_FILENAME' -O '$GENOME_COMP_OUTDIR
    echo $execstr
    $execstr
    RESTRICTION_GENOME_COMP=`expr $RESTRICTION_GENOME_COMP + 1`
fi


# Execution time log file
exectimefile=$OUTPUT_DIR'/ExecTime'
mkdir -p $OUTPUT_DIR
echo 'Execution time for ---> '$DATASET > $exectimefile

#------------------------------------------------
# This script call the Derive_Features_Vector.py
#------------------------------------------------
if [ $MODE_2 == 1 ]
then
  start=`date +%s`
  execstr='./GRAFree.py -d '$DATASET_DIR' -i '$INPUT_FILENAME' -o '$OUTPUT_DIR' -l '$BLOCK_SIZE' -f '$NO_OF_FRAGMENTS' -a '$ASTRAL_PKG' -x '$exectimefile' --redo -w -p -s '$SAVE_DATAFORMAT
  echo $execstr
  $execstr
  end_GRAFree=`date +%s`
  runtime_GRAFree=$((end_GRAFree-start))
  echo '' >> $exectimefile
  echo 'Total Execution time for GRAFree: '$runtime_GRAFree >> $exectimefile
fi


#-----------------------------------------------------------------------------------
# This script compute the bootstrap clade support value wrt the model tree
#-----------------------------------------------------------------------------------
if [ $MODE_3 == 1 ]
then
  echo 'Bootstrap_Support'
  echo '------------------'
  rep=0
  while [ $rep -lt $num_of_replica ]
  do
    BS_OUTPUT_BASEDIR=$OUTPUT_BASE_DIR'/'$DATASET'/Proposed_Bootstrap/'$rep
    BS_DATASET_DIR=$BS_OUTPUT_BASEDIR'/Bootstrapped_Sequences'
    BS_INPUT_FILENAME='Bootstrapped_Replica.fasta'
    BS_OUTPUT_DIR=$OUTPUT_DIR'/Proposed_Bootstrap/'$rep

    regen=true
    if [ "$regen" = true ]
    then
      # Generate replica
      echo '******** #BOOTSTRAP REPLICA: '$rep
      variation_probability=0.3
      execstr='./Proposed_Bootstrapping.py -d '$DATASET_DIR' -o '$BS_DATASET_DIR' -i '$INPUT_FILENAME' -v '$variation_probability
      echo $execstr
      $execstr
    fi
    echo '******** DERIVING BOOTSTRAP TREE: '$rep
    # Derive tree from replica
    execstr='./GRAFree.py -d '$BS_DATASET_DIR' -i '$BS_INPUT_FILENAME' -o '$BS_OUTPUT_DIR' -l '$BLOCK_SIZE' -f '$NO_OF_FRAGMENTS' --redo'
    echo $execstr
    $execstr
    rep=`expr $rep + 1`
  done

fi
