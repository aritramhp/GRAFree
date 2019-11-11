#!/usr/bin/env python

import Header
from Header import *


###-----------------------------------------------------
#''' this function is useful to parse various options for input data processing '''
###-----------------------------------------------------
#def parse_options():  
	#parser = OptionParser()

	#parser.add_option("-I", "--INPTREE", \
		#type="string", \
		#action="store", \
		#dest="INP_TREE", \
		#default="", \
		#help="path of the input tree")	
	
	#parser.add_option("-O", "--OUTPUTTREE", \
		#type="string", \
		#action="store", \
		#dest="OUTPUT_TREE", \
		#default="", \
		#help="path of the output tree")	
		
	#opts, args = parser.parse_args()
	#return opts, args  



#-----------------------------------------------------------------
# This function convert the output of upgma into a newick format
# Output of upgma is a newick file having the info of inner nodes.
#-----------------------------------------------------------------
def ConvertPerfectNewick(input_tree):
	#input_tree = opts.INP_TREE
	#output_tree = opts.OUTPUT_TREE
	
	# tree = dendropy.Tree()
	tree = Tree.get(path=input_tree, schema='newick', preserve_underscores=True)
	tree.write(path=input_tree, schema='newick', unquoted_underscores=True, suppress_internal_node_labels=True,
			   suppress_edge_lengths=True, suppress_rooting=True)
	



#---------------------------------------------------------------------------------------

if __name__=='__main__':
	main()
