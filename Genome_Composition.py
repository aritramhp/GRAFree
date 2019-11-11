#!/usr/bin/env python

import Header
from Header import *

##---------------------------------------------------------------------------------
''' this function is useful to parse various options for input data processing '''
##---------------------------------------------------------------------------------
def parse_options():  
	parser = OptionParser()

	parser.add_option("-D", "--INPDIR", \
		type="string", \
		action="store", \
		dest="INP_DIR", \
		default="", \
		help="path of the input directory containing all species genome in subdirectory")	
	
	parser.add_option("-O", "--OUTDIR", \
		type="string", \
		action="store", \
		dest="OUT_DIR", \
		default="", \
		help="name of the output directory containing all output files in subdirectory (INP_DIR)")	
	
	parser.add_option("-I", "--INPSEQUENCEFILE", \
		type="string", \
		action="store", \
		dest="INP_SEQUENCE_FILE", \
		default="", \
		help="name of the file containing the sequence of the species (all should be same name)")
	
	parser.add_option("-H", "--REVSTRAND", \
		action="store_true", \
		dest="REV_STRAND", \
		default=False, \
		help="If 0, it computes the GFP for forward strands \
					Else, it will computes the GFP for reverse strands")
	
	opts, args = parser.parse_args()
	return opts, args  



#---------------------------------------------------------------------------------------
def main():
	opts, args = parse_options()
	
	INP_DIR = opts.INP_DIR
	OUT_DIR = opts.OUT_DIR
	
	INP_SEQUENCE_FILE = opts.INP_SEQUENCE_FILE
	REV_STRAND = opts.REV_STRAND
	if REV_STRAND == 1:
		REV_STRAND = True
	else:
		REV_STRAND = False
	
	dataset = INP_DIR.split('/')[-1]
	
	if (os.path.isdir(OUT_DIR) == False):
		mkdr_cmd = 'mkdir -p ' + OUT_DIR
		os.system(mkdr_cmd)
	
	OUT_FILE = OUT_DIR + '/' + dataset + '.csv'
	fp_write = open(OUT_FILE,'w')
	fp_write.write('Species Name,Accession Number,Sequence Length,A%,T%,G%,C%,Unrecognized%,AT%,GC%,AT skew, GC skew\n')
	
	rootDir = INP_DIR
		
	for dirName, subdirList, fileList in os.walk(rootDir):
		for fname in fileList:
			if fname == INP_SEQUENCE_FILE:
				input_file = dirName + '/' + INP_SEQUENCE_FILE
				#print'\nFound file: ', input_file
				species_name = dirName.split('/')[-1]
				
				gene_sequences = SeqIO.parse(open(input_file),'fasta')
				A_count = 0
				T_count = 0
				G_count = 0
				C_count = 0
				others = 0
				for gene in gene_sequences:
					seq_id = gene.id
					
					sequence=gene.seq
					
					if REV_STRAND:
						sequence = sequence.reverse_complement()
					
					seq_len = len(sequence)
					
					for i in range(seq_len):
						if sequence[i] == 'A':
							A_count += 1
						elif sequence[i] == 'T':
							T_count += 1
						elif sequence[i] == 'G':
							G_count += 1
						elif sequence[i] == 'C':
							C_count += 1
						else:
							others += 1
					
					A_percent = A_count*100.0/seq_len
					T_percent = T_count*100.0/seq_len
					G_percent = G_count*100.0/seq_len
					C_percent = C_count*100.0/seq_len
					others_percent = others*100.0/seq_len
					
					AT_percent = (A_count+T_count)*100.0/seq_len
					GC_percent = (G_count+C_count)*100.0/seq_len
					
					AT_skew = (A_count-T_count)*1.0/(A_count+T_count)
					GC_skew = (G_count-C_count)*1.0/(G_count+C_count)
					
					fp_write.write(species_name+','+seq_id+','+str(seq_len)+','+str(A_percent)+','+str(T_percent)+','+str(G_percent)+','+str(C_percent)+','+str(others_percent)+','+str(AT_percent)+','+str(GC_percent)+','+str(AT_skew)+','+str(GC_skew)+'\n')
	
	print 'The Genome Composition is written at: ',OUT_FILE
	fp_write.close()	



#---------------------------------------------------------------------------------------

if __name__=='__main__':
	main()
