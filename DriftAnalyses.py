#!/usr/bin/env python

import Header
from Header import *



#----------------------------------
# Compute the distribution of Mean
#----------------------------------
def ComputeMean(x_axis, y_axis):
	xn_array = numpy.array(x_axis)
	yn_array = numpy.array(y_axis)
	xn_mean = []
	yn_mean = []
	
	mean_x = numpy.mean(xn_array)
	mean_y = numpy.mean(yn_array)
	xn_mean.append(mean_x)
	yn_mean.append(mean_y)
	
	return mean_x, mean_y


#---------------------------------
# Compute the Standard Deviation
#---------------------------------
def ComputeSD(x_axis, y_axis):
	xn_array = numpy.array(x_axis)
	yn_array = numpy.array(y_axis)
	xn_std_dev = []
	yn_std_dev = []
	
	std_dev_x = numpy.std(xn_array)
	std_dev_y = numpy.std(yn_array)
	xn_std_dev.append(std_dev_x)
	yn_std_dev.append(std_dev_y)
			
	return std_dev_x, std_dev_y
	
	
#----------------------------------------------------------------
# Compute Eigen values and Eigen vector of the Covariance Matrix
#----------------------------------------------------------------
def ComputeEigenValues(x_axis, y_axis):
	eigen_value = []
	dominant_angle = []
	
	X = numpy.vstack((x_axis, y_axis))
	# Compute the Covariance Matrix
	#print X
	cov_matrix = numpy.cov(X)
	#print cov_matrix
	# Compute eigen value 
	eigval, eigvec = linalg.eig(cov_matrix)
	eigval_list = eigval.tolist()
	eigval = sorted(eigval, key=abs, reverse=True)
	eigen_value.extend(eigval)
	# Compute angle of dominant eigen vector with x-axis
	dominant_pos = eigval_list.index(max(eigval_list, key=abs))
	if eigvec[dominant_pos][0] == 0:
		if eigvec[dominant_pos][1] < 0:
			angle = -90
		else:
			angle = 90
	else:
		angle = math.degrees(math.atan(eigvec[dominant_pos][1] / eigvec[dominant_pos][0]))
	dominant_angle.append(angle)
		
	return eigval, angle


#----------------------------------------------------------------
# Different analyses 
#----------------------------------------------------------------
def DriftAnalyses(x_axis, y_axis, output_dir, analysis_filename, WRITE_INTERMEDIATES):	
	##Max & min
	#x_min = min(x_axis)
	#x_max = max(x_axis)
	#y_min = min(y_axis)
	#y_max = max(y_axis)
	
	# Centroid
	xn_mean, yn_mean = ComputeMean(x_axis, y_axis)
	
	## Standard Deviation
	#xn_std_dev, yn_std_dev = ComputeSD(x_axis, y_axis)
	
	# Compute Eigen values and angle of dominant eigen vector of the Covariance Matrix
	eigenvalues, angle_eigenvector = ComputeEigenValues(x_axis, y_axis)
		
	## Start point and end point
	#mu_S_x = (x_axis[0]-x_min)*1.0/(x_max-x_min)
	#mu_S_y = (y_axis[0]-y_min)*1.0/(y_max-y_min)
	#mu_E_x = (x_axis[-1]-x_min)*1.0/(x_max-x_min)
	#mu_E_y = (y_axis[-1]-y_min)*1.0/(y_max-y_min)
	
	## Horizontal and vertical motion
	#x_temp = sum([abs(x_axis[i]-x_axis[i+1]) for i in range(len(x_axis)-1)])
	#mu_HMN = (x_max - x_min)*1.0/x_temp
	#y_temp = sum([abs(y_axis[i]-y_axis[i+1]) for i in range(len(y_axis)-1)])
	#mu_VMN = (y_max - y_min)*1.0/y_temp
	
	## Aspect ratio
	#mu_AR = (y_max - y_min)*1.0/(math.sqrt((y_max-y_min)**2 + (x_max-x_min)**2))
	
	## Arcedness and straightness
	#arc = len(x_axis) - math.sqrt((y_max-y_min)**2 + (x_max-x_min)**2)
	
	#feat_vect = [x_min, x_max, y_min, y_max, xn_mean, yn_mean, xn_std_dev, yn_std_dev, eigenvalues[0], eigenvalues[1], angle_eigenvector, mu_S_x, mu_S_y, mu_E_x, mu_E_y, mu_HMN, mu_VMN, mu_AR, arc]
	feat_vect = [xn_mean, yn_mean, eigenvalues[0], eigenvalues[1], angle_eigenvector]
		
	#feat_vect = [xn_mean,yn_mean,eigenvalues,angle_eigenvector]
	#print feat_vect
		
	return feat_vect, xn_mean, yn_mean, eigenvalues, angle_eigenvector


#----------------------------------------------------------------------------------------------------
def Analysis(x_cords, y_cords, drift_x, drift_y, NO_OF_FRAGMENTS, output_dir, case, WRITE_INTERMEDIATES):	
	#print '\tFEATURE VECTOR : ', species	
	sub_name = case
	
	#print len(drift_x)
	
	## GFP
	## collect x axis data of gfp	
	#x_axis = x_cords
	## collect y axis data of gfp
	#y_axis = y_cords

	#analysis_filename = 'GFP_Analyses_'+ sub_name
			
	#feat_vect_gfp = DriftAnalyses(x_axis, y_axis, output_dir, analysis_filename, WRITE_INTERMEDIATES)
	#feature_vector.extend(feat_vect_gfp)
		
	# Drift
	drift_feature_vector = []
	drift_xn_mean = []
	drift_yn_mean = []
	drift_eigenvalues = []
	drift_angle_eigenvector = []
	analysis_filename = 'Drift_Analyses_'+ sub_name
	
	for frag in range(NO_OF_FRAGMENTS):
		start = frag * (len(drift_x) / NO_OF_FRAGMENTS)
		end = min((frag+1)*(len(drift_x)/NO_OF_FRAGMENTS), len(drift_x))
		# collect x axis data of Drift	
		frag_x_axis = drift_x[start:end]
		# collect y axis data of Drift
		frag_y_axis = drift_y[start:end]	
				
		feat_vect_drift, xn_mean, yn_mean, eigenvalues, angle_eigenvector = DriftAnalyses(frag_x_axis, frag_y_axis, output_dir, analysis_filename, WRITE_INTERMEDIATES)		
				
		drift_feature_vector.append(feat_vect_drift)
		drift_xn_mean.append(xn_mean)
		drift_yn_mean.append(yn_mean)
		drift_eigenvalues.append(eigenvalues)
		drift_angle_eigenvector.append(angle_eigenvector)
	
	if WRITE_INTERMEDIATES == True:
		output_filename = output_dir + '/' + analysis_filename
		if (os.path.isdir(output_dir) == False):
			mkdr_cmd = 'mkdir -p ' + output_dir
			os.system(mkdr_cmd)
		fp = open(output_filename, 'w')
		#fp.write('X min: ' + str(x_min) + '\n')
		#fp.write('X max: ' + str(x_max) + '\n')
		#fp.write('Y min: ' + str(y_min) + '\n')
		#fp.write('Y max: ' + str(y_max) + '\n')
		fp.write('Mean X: ' + str(drift_xn_mean) + '\n')
		fp.write('Mean Y: ' + str(drift_yn_mean) + '\n')
		#fp.write('Standard Deviation X:  ' + str(xn_std_dev) + '\n')
		#fp.write('Standard Deviation Y:  ' + str(yn_std_dev) + '\n')
		fp.write('Eigenvalues:  ' + str(drift_eigenvalues) + '\n')
		fp.write('Angle eigenvector:  ' + str(drift_angle_eigenvector) + '\n')
		#fp.write('mu_S_x: ' + str(mu_S_x) + '\n')
		#fp.write('mu_S_y: ' + str(mu_S_x) + '\n')		
		#fp.write('mu_E_x: ' + str(mu_E_x) + '\n')
		#fp.write('mu_E_y: ' + str(mu_E_y) + '\n')
		#fp.write('mu_HMN: ' + str(mu_HMN) + '\n')
		#fp.write('mu_VMN: ' + str(mu_VMN) + '\n')
		#fp.write('mu_AR: ' + str(mu_AR) + '\n')
		#fp.write('arc: ' + str(arc) + '\n')		
		fp.close()	
		
	return drift_feature_vector
