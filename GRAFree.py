#!/usr/bin/env python
import Header
from Header import *
import GFP
import Plot
import Drift
import DriftAnalyses
import Compute_Distance_Matrix
import MakeTreeList
import ExecuteASTRAL

##---------------------------------------------------------------------------------
''' this function is useful to parse various options for input data processing '''


##---------------------------------------------------------------------------------
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
    
    parser.add_option("-i", "--INPSEQUENCEFILE",
                      type="string",
                      action="store",
                      dest="INP_SEQUENCE_FILE",
                      default="",
                      help="name of the file containing the sequence of the species (all should be same name)")
    
    parser.add_option("-l", "--BLOCKSIZE",
                      type="int",
                      action="store",
                      dest="BLOCKSIZE",
                      default=0,
                      help="BLOCKSIZE (L)")
    
    parser.add_option("-f", "--FRAGMENT",
                      type="int",
                      action="store",
                      dest="FRAGMENT",
                      default=0,
                      help="FRAGMENT (F)")    
    
    parser.add_option("-a", "--ASTRALPKG",
                      type="string",
                      action="store",
                      dest="ASTRAL_PKG",
                      default="ASTRAL-master-4.10.12/astral.4.10.12.jar",
                      help="Mention the location of the Astral")
    
    parser.add_option("-b", "--BOOTSTRAP",
                      action="store_true",
                      dest="BOOTSTRAP",
                      default=False,
                      help="If TRUE, it execute in BOOTSTRAP mode, where the output directory structutes are different \
                      Else, it will run in normal mode")
    
    parser.add_option("-x", "--TIMELOGFILE",
                      type="string",
                      action="store",
                      dest="TIMELOG_FILE",
                      default="",
                      help="Store the time of executions")
    
    parser.add_option("-w", "--WRITEINTERMEDIATES",
                      action="store_true",
                      dest="WRITE_INTERMEDIATES",
                      default=False,
                      help="If TRUE, it will write the intermediate data in different files")
    
    parser.add_option("-p", "--PLOT",
                      action="store_true",
                      dest="PLOT",
                      default=False,
                      help="If TRUE, it will plot the value in a graph")
    
    parser.add_option("-s", "--SAVEDATAFORMAT",
                      type="string",
                      action="store",
                      dest="SAVE_DATAFORMAT",
                      default="jpg",
                      help="save the plot in this format")
    
    parser.add_option("--redo", "--REDO",
                      action="store_true",
                      dest="REDO",
                      default=False,
                      help="If TRUE, the program executes without taking the intermediate results stored in previous run.")
    
    opts, args = parser.parse_args()
    return opts, args


# ---------------------------------------------------------------------------------------
# Compute GFP
# ---------------------------------------------------------------------------------------
def ComputeGFP(INP_DIR, sp_in_study, INP_SEQUENCE_FILE, CASES, gfp_outdir):
    gfp_file_x = gfp_outdir + '/GFP_x'
    gfp_file_y = gfp_outdir + '/GFP_y'
    GFP_x_dict = dict()
    GFP_y_dict = dict()
    
    # REDO: execute all the results without using the previous results
    if REDO:
        print 'Deriving GFP... '
        # Derive GFP
        for species in sp_in_study:
            inp_seqfile = INP_DIR + '/' + species + '/' + INP_SEQUENCE_FILE
            GFP_x = []
            GFP_y = []
            output_dir = gfp_outdir + '/' + species
            for case in CASES:
                gfp_x, gfp_y = GFP.GraphicalFootPrint(inp_seqfile, output_dir, case, WRITE_INTERMEDIATES, PLOT,
                                                      SAVE_DATAFORMAT)
                GFP_x.append(gfp_x)
                GFP_y.append(gfp_y)
                try:
                    GFP_x_dict.update({species: GFP_x})
                    GFP_y_dict.update({species: GFP_y})
                except:
                    print 'Error in computing GFP. Species: ', species
        # Store dictionary in a file
        if (os.path.isdir(gfp_outdir) == False):
            mkdr_cmd = 'mkdir -p ' + gfp_outdir
            os.system(mkdr_cmd)
        print 'Writing all GFPs in a file...'
        if GFP_x_dict:
            with open(gfp_file_x, 'wb') as dict_items_save:
                cPickle.dump(GFP_x_dict, dict_items_save)
        if GFP_y_dict:
            with open(gfp_file_y, 'wb') as dict_items_save:
                cPickle.dump(GFP_y_dict, dict_items_save)
    else:
        try:
            # Load dictionary from file
            with open(gfp_file_x, 'rb') as dict_items_load:
                print 'Loading GFP file...'
                GFP_x_dict = cPickle.load(dict_items_load)
            with open(gfp_file_y, 'rb') as dict_items_load:
                GFP_y_dict = cPickle.load(dict_items_load)
        except:
            print 'Deriving GFP... '
            # Derive GFP
            for species in sp_in_study:
                inp_seqfile = INP_DIR + '/' + species + '/' + INP_SEQUENCE_FILE
                GFP_x = []
                GFP_y = []
                output_dir = gfp_outdir + '/' + species
                for case in CASES:
                    gfp_x, gfp_y = GFP.GraphicalFootPrint(inp_seqfile, output_dir, case, WRITE_INTERMEDIATES, PLOT,
                                                          SAVE_DATAFORMAT)
                    GFP_x.append(gfp_x)
                    GFP_y.append(gfp_y)
                try:
                    GFP_x_dict.update({species: GFP_x})
                    GFP_y_dict.update({species: GFP_y})
                except:
                    print 'Error in computing GFP. Species: ', species
            # Store dictionary in a file
            if (os.path.isdir(gfp_outdir) == False):
                mkdr_cmd = 'mkdir -p ' + gfp_outdir
                os.system(mkdr_cmd)
            print 'Writing all GFPs in a file...'
            if GFP_x_dict:
                with open(gfp_file_x, 'wb') as dict_items_save:
                    cPickle.dump(GFP_x_dict, dict_items_save)
            if GFP_y_dict:
                with open(gfp_file_y, 'wb') as dict_items_save:
                    cPickle.dump(GFP_y_dict, dict_items_save)
    
    return GFP_x_dict, GFP_y_dict


# ---------------------------------------------------------------------------------------
# Compute drift for the sequence 
# ---------------------------------------------------------------------------------------
def ComputeDrift(sp_in_study, GFP_x_dict, GFP_y_dict, drift_outdir, BLOCKSIZE, CASES):
    drift_file_x = drift_outdir + '/Drift_x'
    drift_file_y = drift_outdir + '/Drift_y'
    Drift_x_dict = dict()
    Drift_y_dict = dict()
    
    # REDO: execute all the results without using the previous results
    if REDO:
        print 'Deriving Drift... '
        # Derive Drift
        for species in sp_in_study:
            Drift_x = []
            Drift_y = []
            output_dir = drift_outdir + '/' + species
            for c in range(len(CASES)):
                case = CASES[c]
                GFP_x = GFP_x_dict[species][c]
                GFP_y = GFP_y_dict[species][c]
                dft_x, dft_y = Drift.Drift(GFP_x, GFP_y, output_dir, BLOCKSIZE, case, WRITE_INTERMEDIATES, PLOT,
                                           SAVE_DATAFORMAT)
                Drift_x.append(dft_x)
                Drift_y.append(dft_y)
            try:
                Drift_x_dict.update({species: Drift_x})
                Drift_y_dict.update({species: Drift_y})
            except:
                print 'Error in computing Drift. Species: ', species
        # Store dictionary in a file
        if os.path.isdir(drift_outdir) == False:
            mkdr_cmd = 'mkdir -p ' + drift_outdir
            os.system(mkdr_cmd)
        print 'Writing all Drifts in a file...'
        if Drift_x_dict:
            with open(drift_file_x, 'wb') as dict_items_save:
                cPickle.dump(Drift_x_dict, dict_items_save)
        if Drift_y_dict:
            with open(drift_file_y, 'wb') as dict_items_save:
                cPickle.dump(Drift_y_dict, dict_items_save)
    else:
        try:
            # Load dictionary from file
            with open(drift_file_x, 'rb') as dict_items_load:
                print 'Loading Drift file...'
                Drift_x_dict = cPickle.load(dict_items_load)
            with open(drift_file_y, 'rb') as dict_items_load:
                Drift_y_dict = cPickle.load(dict_items_load)
        except:
            print 'Deriving Drift... '
            # Derive Drift
            for species in sp_in_study:
                Drift_x = []
                Drift_y = []
                output_dir = drift_outdir + '/' + species
                for c in range(len(CASES)):
                    case = CASES[c]
                    GFP_x = GFP_x_dict[species][c]
                    GFP_y = GFP_y_dict[species][c]
                    dft_x, dft_y = Drift.Drift(GFP_x, GFP_y, output_dir, BLOCKSIZE, case, WRITE_INTERMEDIATES, PLOT,
                                               SAVE_DATAFORMAT)
                    Drift_x.append(dft_x)
                    Drift_y.append(dft_y)
                try:
                    Drift_x_dict.update({species: Drift_x})
                    Drift_y_dict.update({species: Drift_y})
                except:
                    print 'Error in computing Drift. Species: ', species
            # Store dictionary in a file
            if os.path.isdir(drift_outdir) == False:
                mkdr_cmd = 'mkdir -p ' + drift_outdir
                os.system(mkdr_cmd)
            print 'Writing all Drifts in a file...'
            if Drift_x_dict:
                with open(drift_file_x, 'wb') as dict_items_save:
                    cPickle.dump(Drift_x_dict, dict_items_save)
            if Drift_y_dict:
                with open(drift_file_y, 'wb') as dict_items_save:
                    cPickle.dump(Drift_y_dict, dict_items_save)
    
    return Drift_x_dict, Drift_y_dict


# ---------------------------------------------------------------------------------------
# Compute FEATURE_VECTOR
# ---------------------------------------------------------------------------------------
def ComputeFV(sp_in_study, Drift_x_dict, Drift_y_dict, FRAGMENT, feat_outdir, CASES):
    featvect_file = feat_outdir + '/FeatVect'
    FeatVect_dict = dict()
    
    # REDO: execute all the results without using the previous results
    if REDO:
        print 'Deriving Feature Vector... '
        # Derive Drift
        for species in sp_in_study:
            feat_vect = []
            output_dir = feat_outdir + '/' + species
            for c in range(len(CASES)):
                case = CASES[c]
                GFP_x = GFP_y = []
                DRIFT_x = Drift_x_dict[species][c]
                DRIFT_y = Drift_y_dict[species][c]
                fv = DriftAnalyses.Analysis(GFP_x, GFP_y, DRIFT_x, DRIFT_y, FRAGMENT, output_dir, case,
                                            WRITE_INTERMEDIATES)
                feat_vect.append(fv)
            try:
                FeatVect_dict.update({species: feat_vect})
            except:
                print 'Error in computing Feature Vector. Species: ', species
        # Store dictionary in a file
        if os.path.isdir(feat_outdir) == False:
            mkdr_cmd = 'mkdir -p ' + feat_outdir
            os.system(mkdr_cmd)
        print 'Writing all Feature Vector in a file...\n'
        if FeatVect_dict:
            with open(featvect_file, 'wb') as dict_items_save:
                cPickle.dump(FeatVect_dict, dict_items_save)
    else:
        try:
            # Load dictionary from file
            with open(featvect_file, 'rb') as dict_items_load:
                print 'Loading Features file...'
                FeatVect_dict = cPickle.load(dict_items_load)
        except:
            print 'Deriving Feature Vector... '
            # Derive Drift
            for species in sp_in_study:
                feat_vect = []
                output_dir = feat_outdir + '/' + species
                for c in range(len(CASES)):
                    case = CASES[c]
                    GFP_x = GFP_y = []
                    DRIFT_x = Drift_x_dict[species][c]
                    DRIFT_y = Drift_y_dict[species][c]
                    fv = DriftAnalyses.Analysis(GFP_x, GFP_y, DRIFT_x, DRIFT_y, FRAGMENT, output_dir, case,
                                                WRITE_INTERMEDIATES)
                    feat_vect.append(fv)
                try:
                    FeatVect_dict.update({species: feat_vect})
                except:
                    print 'Error in computing Feature Vector. Species: ', species
            # Store dictionary in a file
            if os.path.isdir(feat_outdir) == False:
                mkdr_cmd = 'mkdir -p ' + feat_outdir
                os.system(mkdr_cmd)
            print 'Writing all Feature Vector in a file...\n'
            if FeatVect_dict:
                with open(featvect_file, 'wb') as dict_items_save:
                    cPickle.dump(FeatVect_dict, dict_items_save)
    
    return FeatVect_dict


# ------------------------------------------------------------------
# Compute the Distance matrix by calling Compute_Distance_Matrix.py
# ------------------------------------------------------------------
def DISTANCE_ProcessCall(FEATURE_VECTOR, SPECIES_NAME, BLOCKSIZE, FRAGMENT, alpha, case, DIST_FUNC, output_dir):
    if DIST_FUNC == 'GRAFree':
        distance_matrix = Compute_Distance_Matrix.ComputeDistance(SPECIES_NAME, FEATURE_VECTOR, alpha, FRAGMENT,
                                                                  output_dir, case, DIST_FUNC, WRITE_INTERMEDIATES)
    elif DIST_FUNC == 'Euclidean':
        distance_matrix = Compute_Distance_Matrix.ComputeDistance(SPECIES_NAME, FEATURE_VECTOR, alpha, FRAGMENT,
                                                                  output_dir, case, DIST_FUNC, WRITE_INTERMEDIATES)
    return distance_matrix


# ---------------------------------------------------------------------------------------
def main():
	opts, args = parse_options()

	INP_DIR = opts.INP_DIR
	global OUT_DIR
	OUT_DIR = opts.OUT_DIR
	INP_SEQUENCE_FILE = opts.INP_SEQUENCE_FILE

	CASES = ['1', '2', '3']
	BLOCKSIZE = opts.BLOCKSIZE
	FRAGMENT = opts.FRAGMENT
	ALPHA_LIST = [0.50]
	if BLOCKSIZE <= 0 or FRAGMENT <= 0:
		print 'Blocksize and number of fragments must be a positive integer.'
		sys.exit()

	ASTRAL_PKG = opts.ASTRAL_PKG

	TIMELOG_FILE = opts.TIMELOG_FILE

	global WRITE_INTERMEDIATES
	WRITE_INTERMEDIATES = opts.WRITE_INTERMEDIATES
	global PLOT
	PLOT = opts.PLOT
	global SAVE_DATAFORMAT
	SAVE_DATAFORMAT = opts.SAVE_DATAFORMAT

	global DATASET
	DATASET = INP_DIR.split('/')[-1]
	global REDO
	REDO = opts.REDO

	BOOTSTRAP = opts.BOOTSTRAP
	if BOOTSTRAP == True:
		log = True
	else:
		log = False

	if log == True or not TIMELOG_FILE == '':
		fp_time = open(TIMELOG_FILE, 'a')

	list_of_sp = os.listdir(INP_DIR)

	# Compute GFP by calling the function ComputeGFP
	print '::::: Computing GFPs :::::'
	gfp_start = time.time()
	gfp_outdir = OUT_DIR + '/GFP'
	GFP_x_dict, GFP_y_dict = ComputeGFP(INP_DIR, list_of_sp, INP_SEQUENCE_FILE, CASES, gfp_outdir)
	gfp_end = time.time()
	gfp_exectime = gfp_end - gfp_start
	if log == True or not TIMELOG_FILE == '':
		fp_time.write('GFP : ' + str(gfp_exectime) + '\n')

	# Compute Drift by calling the function ComputeDrift
	print '::::: Computing Drifts :::::'
	# Derive Drift
	drift_outdir = OUT_DIR + '/Drift/L_' + str(BLOCKSIZE)
	drift_start = time.time()
	Drift_x_dict, Drift_y_dict = ComputeDrift(list_of_sp, GFP_x_dict, GFP_y_dict, drift_outdir, BLOCKSIZE, CASES)
	drift_end = time.time()
	drift_exectime = drift_end - drift_start
	if log == True or not TIMELOG_FILE == '':
		fp_time.write('\nDRIFT: L=' + str(BLOCKSIZE) + ' : ' + str(drift_exectime) + '\n')

	# Compute FEAT_VECT by calling the function ComputeFV
	print '::::: Deriving Feature Vectors :::::'
	# Derive feature vector
	feat_outdir = OUT_DIR + '/FeatVect/L_' + str(BLOCKSIZE) + '/F_' + str(FRAGMENT)
	feat_start = time.time()
	FeatVect_dict = ComputeFV(list_of_sp, Drift_x_dict, Drift_y_dict, FRAGMENT, feat_outdir, CASES)
	feat_end = time.time()
	feat_exectime = feat_end - feat_start
	if log == True or not TIMELOG_FILE == '':
		fp_time.write(
			'\tFeature Vector: L=' + str(BLOCKSIZE) + ', F=' + str(FRAGMENT) + ' : ' + str(feat_exectime) + '\n')

	# Compute DISTANCE MATRIX by calling the function DISTANCE_ProcessCall
	FEATURE_VECTOR = []
	for cs in range(len(CASES)):
		feat = []
		for sp in list_of_sp:
			feat.append(FeatVect_dict[sp][cs])
		FEATURE_VECTOR.append(feat)

	print '::::: Deriving Distance Matrices and UPGMA :::::'
	for cs in range(len(CASES)):
		case = CASES[cs]
		if log == True or not TIMELOG_FILE == '':
			fp_time.write('Case: ' + case + ' :\n')        
		for alpha in ALPHA_LIST:
			alpha_dir = str(int(alpha * 100))
			dm_start = time.time()
			output_dir = OUT_DIR + '/OUTPUT_TREE/GRAFree/L_' + str(BLOCKSIZE) + '/F_' + str(
				FRAGMENT) + '/' + alpha_dir
			#distace_matrix = DISTANCE_ProcessCall(FEATURE_VECTOR[cs], list_of_sp, BLOCKSIZE, FRAGMENT, alpha,
													#case, DISTANCE_FUNCTION, output_dir)
			distance_matrix = Compute_Distance_Matrix.ComputeDistance(list_of_sp, FEATURE_VECTOR[cs], alpha, FRAGMENT,
                                                                  output_dir, case, WRITE_INTERMEDIATES)
			dm_end = time.time()
			dm_exectime = dm_end - dm_start
			if log == True or not TIMELOG_FILE == '':
				fp_time.write(
					'\tL=' + str(BLOCKSIZE) + ', F=' + str(FRAGMENT) + ', A=' + str(alpha) + ' : ' + str(
						dm_exectime) + '\n')
			# derive tree using UPGMA
			upgma_start = time.time()
			Compute_Distance_Matrix.CallUPGMA(list_of_sp, distance_matrix, output_dir, case)
			upgma_end = time.time()
			upgma_exectime = upgma_end - upgma_start
			if log == True or not TIMELOG_FILE == '':
				fp_time.write(
					'\tUPGMA: L=' + str(BLOCKSIZE) + ', F=' + str(FRAGMENT) + ', A=' + str(alpha) + ': ' +
					str(upgma_exectime) + '\n')        
	FEATURE_VECTOR = None
	gc.collect()
	# --------------------------------------------------------------------------------------------------------------------
	# Merge tree
	# --------------------------------------------------------------------------------------------------------------------
	for alpha in ALPHA_LIST:
		alpha_dir = str(int(alpha * 100))
		# Make treelist
		INP_DIR = OUT_DIR + '/OUTPUT_TREE/GRAFree/L_' + str(BLOCKSIZE) + '/F_' + str(FRAGMENT) + '/' + alpha_dir
		MakeTreeList.TreeListing(INP_DIR, CASES)
		# Call ASTRAL
		print '::::: Executing ASTRAL :::::'
		astral_start = time.time()
		ExecuteASTRAL.main(INP_DIR, ASTRAL_PKG)
		astral_end = time.time()
		astral_exectime = astral_end - astral_start
		if log == True or not TIMELOG_FILE == '':
			fp_time.write('\n\t***ASTRAL: L=' + str(BLOCKSIZE) + ', F=' + str(FRAGMENT) + ', A=' + str(
				alpha) + ' : ' + str(astral_exectime) + '\n')
	DRIFT_coord_x = None
	gc.collect()
	DRIFT_coord_y = None
	gc.collect()

	if log == True or not TIMELOG_FILE == '':
		fp_time.close()


# ---------------------------------------------------------------------------------------

if __name__ == '__main__':
    main()
