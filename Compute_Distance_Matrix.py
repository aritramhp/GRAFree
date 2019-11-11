#!/usr/bin/env python

import Header
from Header import *
import TreeConstruction


# -----------------------------------------------------------------
# This function convert the output of upgma into a newick format
# Output of upgma is a newick file having the info of inner nodes.
# -----------------------------------------------------------------
def ConvertPerfectNewick(UPGMA_output, output_dir, sub_name):
	tree = Tree.get(path=UPGMA_output, schema='newick', preserve_underscores=True)
	treefile_name = output_dir + '/tree_' + sub_name + '.tre'
	treefile_name_weighted = output_dir + '/tree_' + sub_name + '_weighted.tre'
	tree.write(path=treefile_name, schema='newick', unquoted_underscores=True, suppress_internal_node_labels=True,
				suppress_edge_lengths=True,
				suppress_rooting=False)  # write tree w/o edge lengths
	tree.write_to_path(treefile_name_weighted, 'newick', suppress_internal_node_labels=True,
						suppress_rooting=True)  # write tree with edge lengths


# ---------------------------------------------------------------
# Call UPGMA
# ---------------------------------------------------------------
def CallUPGMA(names, dist_matrix, output_dir, sub_name):
	m = TreeConstruction._DistanceMatrix(names, dist_matrix)
	constructor = TreeConstruction.DistanceTreeConstructor()
	tree = constructor.upgma(m)
	UPGMA_output = output_dir + "/TreeUPGMA.newick"
	Phylo.write(tree, UPGMA_output, 'newick')
	print "Tree saved as TreeUPGMA.newick in: ", output_dir
	ConvertPerfectNewick(UPGMA_output, output_dir, sub_name)


# ---------------------------------------------------------------------------
# This function compute the distance matrix using the following formula
# (proposed by JM):
# for each fragment, dist = alpha * sqrt(mean1^2 + mean2^2 - 2*mean1*mean2*cos(theta1-theta2) +
# 			(1-alpha) * sqrt((eig1_1 - eig2_1)^2 + (eig1_2 - eig2_2)^2)
# Total dist = sum of all dist.
# ---------------------------------------------------------------------------
def ComputeDistanceMatrix_GRAFree(FEAT_VECT, ALPHA, FRAGMENTS):
	no_of_sp = len(FEAT_VECT)

	DISTANCE_MATRIX = []
	for i in range(no_of_sp):
		DISTANCE_MATRIX.append([])
		
		for j in range(i):
			frag_dist = []
			for frag in range(FRAGMENTS):
				mean_x_i = FEAT_VECT[i][frag][0]
				mean_y_i = FEAT_VECT[i][frag][1]
				eig_root_i = [FEAT_VECT[i][frag][2], FEAT_VECT[i][frag][3]]
				theta_i = FEAT_VECT[i][frag][4]
				
				mean_x_j = FEAT_VECT[j][frag][0]
				mean_y_j = FEAT_VECT[j][frag][1]
				eig_root_j = [FEAT_VECT[j][frag][2], FEAT_VECT[j][frag][3]]
				theta_j = FEAT_VECT[j][frag][4]
				
				mean_1_sqr = numpy.dot([mean_x_i, mean_y_i], [mean_x_i, mean_y_i])
				mean_2_sqr = numpy.dot([mean_x_j, mean_y_j], [mean_x_j, mean_y_j])
				theta_diff = abs(theta_i - theta_j)
				d1 = mean_1_sqr + mean_2_sqr
				d2 = 2 * numpy.dot([mean_x_i, mean_y_i], [mean_x_j, mean_y_j]) * (math.cos(math.radians(theta_diff)))
				d3 = math.sqrt(pow((eig_root_i[0] - eig_root_j[0]), 2) + pow((eig_root_i[1] - eig_root_j[1]), 2))
				d4 = math.sqrt(d1 - d2)
				dist = ALPHA * d4 + (1 - ALPHA) * d3
				frag_dist.append(dist)
			
			distn_array = numpy.array(frag_dist)
			distance = numpy.mean(distn_array)
			
			DISTANCE_MATRIX[-1].append(distance)
		
		DISTANCE_MATRIX[-1].append(0)
	# print DISTANCE_MATRIX
	return DISTANCE_MATRIX


# ---------------------------------------------------------------------------
# This function compute the distance matrix using the euclidean distance
# for each fragment, dist = sqrt((mean_x_1 - mean_x_2)^2 + (mean_y_1 - mean_y_2)^2
# 						+ (eig_1 - eig_2)^2 + (theta1 - theta2)^2)
# Total dist = sum of all dist.
# ---------------------------------------------------------------------------
def ComputeDistanceMatrix_Euclidean(FEAT_VECT, ALPHA, FRAGMENTS):
	no_of_sp = len(FEAT_VECT)

	# print zero_feature
	DISTANCE_MATRIX = []
	for i in range(no_of_sp):
		DISTANCE_MATRIX.append([])
		for j in range(i):
			frag_dist = []
			for frag in range(FRAGMENTS):
				mean_x_i = FEAT_VECT[i][frag][0]
				mean_y_i = FEAT_VECT[i][frag][1]
				eig_root_i = [FEAT_VECT[i][frag][2], FEAT_VECT[i][frag][3]]
				theta_i = FEAT_VECT[i][frag][4]
				
				mean_x_j = FEAT_VECT[j][frag][0]
				mean_y_j = FEAT_VECT[j][frag][1]
				eig_root_j = [FEAT_VECT[j][frag][2], FEAT_VECT[j][frag][3]]
				theta_j = FEAT_VECT[j][frag][4]
				
				dist = (mean_x_i - mean_x_j) ** 2 + (mean_y_i - mean_y_j) ** 2 + (
						eig_root_i[0] - eig_root_j[0]) ** 2 + (eig_root_i[1] - eig_root_j[1]) ** 2 + (
								theta_i - theta_j) ** 2
				dist = math.sqrt(dist)
				frag_dist.append(dist)
			
			distn_array = numpy.array(frag_dist)
			distance = numpy.mean(distn_array)
			DISTANCE_MATRIX[-1].append(distance)
		
		DISTANCE_MATRIX[-1].append(0)

	return DISTANCE_MATRIX


# ---------------------------------------------------------------
def ComputeDistance(SPECIES_NAME, FEAT_VECT, ALPHA, FRAGMENT, output_dirname, sub_name, WRITE_INTERMEDIATES):
	no_of_sp = len(SPECIES_NAME)

	DISTANCE_MATRIX = ComputeDistanceMatrix_GRAFree(FEAT_VECT, ALPHA, FRAGMENT)

	if os.path.isdir(output_dirname) == False:
		mkdr_cmd = 'mkdir -p ' + output_dirname
		os.system(mkdr_cmd)
    
	# write the distance matrix in a file (.csv)
	if WRITE_INTERMEDIATES == True:
		output_Dist_raw = output_dirname + '/Raw_Distance_matrix_' + sub_name
		fp = open(output_Dist_raw, 'w')
		fp.write(str(DISTANCE_MATRIX))
		fp.close()
		
		output_Dist = output_dirname + '/Distance_matrix_' + sub_name
		fp = open(output_Dist, 'w')
		for i in range(no_of_sp):
			fp.write(SPECIES_NAME[i] + ':\t\t')
			fp.write(str(DISTANCE_MATRIX[i]) + '\n')
		fp.close()
		
		output_Sp = output_dirname + '/Sp_name_' + sub_name
		fp = open(output_Sp, 'w')
		fp.write(str(SPECIES_NAME))
		fp.close()

	return DISTANCE_MATRIX
