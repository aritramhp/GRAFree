#!/usr/bin/env python

import Header
from Header import *
import Plot



#---------------------------------------------------------------------------------------
# Compute GRAPHICAL-FOOT-PRINT OF DNA
#---------------------------------------------------------------------------------------
def GraphicalFootPrint(filename, output_dir, case, WRITE_INTERMEDIATES, PLOT, SAVE_DATAFORMAT):
	#Compute the footprint of the gene sequence
	gene_sequences = SeqIO.parse(open(filename),'fasta')
	for gene in gene_sequences:
		sequence=gene.seq				
		seq_len = len(sequence)
		gfp_xn = []
		gfp_yn = []
		xn=[0]
		yn=[0]
		
		# CASE--1 (Nandy '94)
		if case == '1':
			for i in range(seq_len):
				if sequence[i] == 'G':
					xn.append(xn[-1]+1)
					yn.append(yn[-1])
				elif sequence[i] == 'A':
					xn.append(xn[-1]-1)
					yn.append(yn[-1])
				elif sequence[i] == 'C':
					xn.append(xn[-1])
					yn.append(yn[-1]+1)
				elif sequence[i] == 'T':
					xn.append(xn[-1])
					yn.append(yn[-1]-1)
				else:
					xn.append(xn[-1])
					yn.append(yn[-1])
			sub_name = '_1'
		#CASE--2 (Gates '86)
		elif case == '2':
			for i in range(seq_len):
				if sequence[i] == 'C':
					xn.append(xn[-1]+1)
					yn.append(yn[-1])
				elif sequence[i] == 'G':
					xn.append(xn[-1]-1)
					yn.append(yn[-1])
				elif sequence[i] == 'T':
					xn.append(xn[-1])
					yn.append(yn[-1]+1)
				elif sequence[i] == 'A':
					xn.append(xn[-1])
					yn.append(yn[-1]-1)
				else:
					xn.append(xn[-1])
					yn.append(yn[-1])
			sub_name = '_2'
		# CASE--3 (Leong & Morgenthaler '95)
		elif case == '3':
			for i in range(seq_len):
				if sequence[i] == 'A':
					xn.append(xn[-1]+1)
					yn.append(yn[-1])
				elif sequence[i] == 'C':
					xn.append(xn[-1]-1)
					yn.append(yn[-1])
				elif sequence[i] == 'T':
					xn.append(xn[-1])
					yn.append(yn[-1]+1)
				elif sequence[i] == 'G':
					xn.append(xn[-1])
					yn.append(yn[-1]-1)
				else:
					xn.append(xn[-1])
					yn.append(yn[-1])
			sub_name = '_3'
			
		gfp_xn.append(xn)
		gfp_yn.append(yn)			
		
		# write the locus of the footprint in a file (.csv)
		if WRITE_INTERMEDIATES == True:				
			if (os.path.isdir(output_dir) == False):
				mkdr_cmd = 'mkdir -p ' + output_dir
				os.system(mkdr_cmd)
			output_FP_X=output_dir+'/xn'+sub_name
			fp=open(output_FP_X,'w')
			fp.write(str(xn))
			fp.close()
			output_FP_Y=output_dir+'/yn'+sub_name
			fp=open(output_FP_Y,'w')
			fp.write(str(yn))
			fp.close()
			print '\tWrite file : ', output_FP_X
			
		# Plot the GFP by calling Plot.py
		if PLOT == True:
			output_plotfile = output_dir + '/GFP_' + case + '.' + SAVE_DATAFORMAT
			Plot.mainPlotting(xn, yn, output_plotfile)
			print '\tPlotted the GFP at: ', output_plotfile
	return xn, yn