#!/usr/bin/env python

import Header
from Header import *
import Plot

#----------------------------
# Set debug level
#----------------------------
#DEBUG_LEVEL = 3


#-----------------------------------------------------------------
# Compute the Drift - the differece between start and end points 
# of the window
#-----------------------------------------------------------------
def ComputeDrift(x_axis, y_axis):
	#print x_axis[-1], x_axis[0]
	drift_x = x_axis[-1] - x_axis[0]
	drift_y = y_axis[-1] - y_axis[0]
	return drift_x, drift_y



def Drift(GFP_x, GFP_y, output_dir, block_size, case, WRITE_INTERMEDIATES, PLOT, SAVE_DATAFORMAT):
	drift_x = []
	drift_y = []
	
	case = int(case)
	
	# collect x axis data
	x_axis = GFP_x
	# collect y axis data
	y_axis = GFP_y

	dft_x = []
	dft_y = []

	for w in range(len(x_axis)-block_size+1):
		window =  min((len(x_axis) - w), block_size)
		d_x, d_y = ComputeDrift(x_axis[w : w+window], y_axis[w : w+window])		
		dft_x.append(d_x)
		dft_y.append(d_y)

	drift_x.append(dft_x)
	drift_y.append(dft_y)
	
	# write the locus of the footprint in a file (.csv)
	if WRITE_INTERMEDIATES == True:
		output_dirname = output_dir
		if (os.path.isdir(output_dirname) == False):
			mkdr_cmd = 'mkdir -p ' + output_dirname
			os.system(mkdr_cmd)
		output_FP_X = output_dirname+'/Drift_xn'+str(case)
		fp=open(output_FP_X,'w')
		fp.write(str(dft_x))
		fp.close()
		print '\tDrift have stored at: ', output_FP_X
		output_FP_Y = output_dirname+'/Drift_yn'+str(case)
		fp=open(output_FP_Y,'w')
		fp.write(str(dft_y))
		fp.close()
		
		print '\tDrift have stored at: ', output_FP_Y
	
	# Plot the GFP by calling Plot.py
		if PLOT == True:
			output_plotfile = output_dir + '/Drift_' + str(case) + '.' + SAVE_DATAFORMAT
			Plot.mainPlotting(dft_x, dft_y, output_plotfile, block_size)
			print '\tPlotted the Drift at: ', output_plotfile
	
	return drift_x[0], drift_y[0]
