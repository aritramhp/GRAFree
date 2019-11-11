#!/usr/bin/env python

import Header
from Header import *


def parse_options():
    parser = OptionParser()
    
    parser.add_option("-d", "--INPDIR",
                      type="string",
                      action="store",
                      dest="INP_DIR",
                      default="",
                      help="path of the input directory containing all species genome in subdirectory")
    
    parser.add_option("-o", "--OUTDIR",
                      type="string",
                      action="store",
                      dest="OUT_DIR",
                      default="",
                      help="name of the output directory containing all output files in subdirectory (INP_DIR)")
    
    parser.add_option("-i", "--INPFILENAME",
                      type="string",
                      action="store",
                      dest="INP_FILENAME",
                      default="",
                      help="name of the file containing the gene sequences of the species (all should be same name)")
    
    parser.add_option("-v", "--VARIATION",
                      type="float",
                      action="store",
                      dest="variation_probability",
                      default="0.3",
                      help="probability of the variation")
    
    opts, args = parser.parse_args()
    return opts, args


# ------------------------------------------------------------------------------------
# Apply mutation with a fixed uniform probability
# ------------------------------------------------------------------------------------
def Mutation(sequence, variation_probability):
    sequence = list(sequence)
    probability_scale = variation_probability * 10.0 / 100
    problimit_mutation = probability_scale / 3
    problimit_insert = probability_scale * 2 / 3
    problimit_delete = probability_scale * 3 / 3
    
    modified_sequence = []
    
    ch = 0
    while ch < len(sequence):
        residues = ['A', 'T', 'G', 'C']
        
        random_change = random.uniform(0, 10)  # probability to mutate a particular residue
        # Mutation
        if random_change <= problimit_mutation:
            if sequence[ch] in residues:
                random_select = random.uniform(0, 3)  # probability to change the residue
                residues.remove(sequence[ch])
                if random_select >= 0 and random_select <= 1:
                    modified_sequence.append(residues[0])
                elif random_select > 1 and random_select <= 2:
                    modified_sequence.append(residues[1])
                elif random_select > 2 and random_select <= 3:
                    modified_sequence.append(residues[2])
            ch += 1
        # Insertion
        elif random_change > problimit_mutation and random_change <= problimit_insert:
            random_select = random.uniform(0, 4)
            if random_select >= 0 and random_select <= 1:
                new = residues[0]
            elif random_select > 1 and random_select <= 2:
                new = residues[1]
            elif random_select > 2 and random_select <= 3:
                new = residues[2]
            elif random_select > 3 and random_select <= 4:
                new = residues[3]
            modified_sequence.append(new)
        # Deletion
        elif random_change > problimit_insert and random_change <= problimit_delete:
            ch += 1
        # No change
        else:
            modified_sequence.append(sequence[ch])
            ch += 1
    
    modified_sequence = ''.join(modified_sequence)
    return modified_sequence


# ------------------------------------------------------------------------------------
# Bootstrapping
# ------------------------------------------------------------------------------------
def Bootstrap(input_file, output_dir, variation_probability):
    # variation_probability = 0.3		# % of variation within intraspecific genome
    fp_input_file = open(input_file, 'r')
    first_line = fp_input_file.readline()[1:]
    record_id = first_line
    
    lines = fp_input_file.readlines()
    record_seq = ''.join(lines)
    record_seq = record_seq.replace('\n', '')
    
    mutated_seq = Mutation(record_seq, variation_probability)
    record = SeqRecord(Seq(mutated_seq), record_id, '', '')
    
    # Write the replica in a fasta file
    species_name = input_file.split('/')[-2]
    replica_dir = output_dir + '/' + species_name
    if (os.path.isdir(replica_dir) == False):
        mkdr_cmd = 'mkdir -p ' + replica_dir
        os.system(mkdr_cmd)
    replica_file = replica_dir + '/Bootstrapped_Replica.fasta'
    SeqIO.write(record, replica_file, "fasta")


# ---------------------------------------------------------------------------------------
# This function bootstrap a single sequence according to our proposed method
# ---------------------------------------------------------------------------------------
def main():
    opts, args = parse_options()
    
    input_dir = opts.INP_DIR
    output_dir = opts.OUT_DIR
    inputfilename = opts.INP_FILENAME
    variation_probability = opts.variation_probability
    
    if (os.path.isdir(output_dir) == False):
        mkdr_cmd = 'mkdir -p ' + output_dir
        os.system(mkdr_cmd)
    
    rootDir = input_dir
    for dirName, subdirList, fileList in os.walk(rootDir):
        for fname in fileList:
            if fname == inputfilename:
                input_file = dirName + '/' + inputfilename
                print "input_file: ", input_file
                Bootstrap(input_file, output_dir, variation_probability)


# ---------------------------------------------------------------------------------------

if __name__ == '__main__':
    main()
