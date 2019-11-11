#!/usr/bin/env python

import Header
from Header import *
import MakePerfectNewickTree



##-----------------------------------------------------
''' this function is useful to parse various options for input data processing '''
##-----------------------------------------------------
def parse_options():  
    parser = OptionParser()

    parser.add_option("-D", "--INPDIR", \
        type="string", \
        action="store", \
        dest="INP_DIR", \
        default="", \
        help="path of the input directory containing all species genome in subdirectory")
    
    parser.add_option("-R", "--PCKGBASEDIR", \
        type="string", \
        action="store", \
        dest="PCKG_BASE_DIR", \
        default="", \
        help="root path of the packages used here (here root location of COSPED)")
        
    opts, args = parser.parse_args()
    return opts, args


# -----------------------------------------------------------------
# This function convert the output of upgma into a newick format
# Output of upgma is a newick file having the info of inner nodes.
# -----------------------------------------------------------------
def ConvertPerfectNewick(astraltreefile):
    tree = Tree.get(path=astraltreefile, schema='newick', preserve_underscores=True)
    # treefile_name_weighted = astraltreefile + 'weighted.tre'
    tree.write(path=astraltreefile, schema='newick', unquoted_underscores=True, suppress_internal_node_labels=True,
               suppress_edge_lengths=True,
               suppress_rooting=False)  # write tree w/o edge lengths
    # tree.write_to_path(treefile_name_weighted, 'newick', suppress_internal_node_labels=True,
    #                    suppress_rooting=True)  # write tree with edge lengths
    
    
#----------------------------------------------------------------------------------------------
# This program execute ASTRAL for all TreeList.newick file found in the specific location
#----------------------------------------------------------------------------------------------
def main(input_dir,astralexec):
    rootDir = input_dir
    #print 'Start executing ASTRAL-2: ', input_dir
    
    treelistfile = 'TreeList.newick'
    print 'Start executing ASTRAL-2: ', input_dir
    filename = rootDir + '/' + treelistfile
    outputtreedir=rootDir + '/ASTRAL2'
    
    if (os.path.isfile(filename) == True):
        if (os.path.isdir(outputtreedir) == False):
            mkdr_cmd = 'mkdir -p ' + outputtreedir
            os.system(mkdr_cmd)
        
        astraltreefile = outputtreedir + '/ASTRAL_newick.tre'
        outtextfile = outputtreedir + '/Complete_Desription.txt'
        # Call ASTRAL
        astralcmd = 'java -jar ' + astralexec + ' -i ' + filename + ' -o ' + astraltreefile + ' > ' + outtextfile
        os.system(astralcmd)

        ConvertPerfectNewick(astraltreefile)
    else:
        print 'File not Found: ', filename
    
    
    
    #for dirName, subdirList, fileList in os.walk(rootDir):
        #for fname in fileList:
            #if fname == treelistfile:
                #print 'Start executing ASTRAL-2: ', input_dir
                #filename = dirName + '/' + fname
                #outputtreedir=dirName + '/ASTRAL2'
                #if (os.path.isdir(outputtreedir) == False):
                    #mkdr_cmd = 'mkdir -p ' + outputtreedir
                    #os.system(mkdr_cmd)
                
                #astraltreefile = outputtreedir + '/ASTRAL_newick.tre'
                #outtextfile = outputtreedir + '/Complete_Desription.txt'
                ## Call ASTRAL
                #astralexec=pkg_basedir + '/Species_Tree_Estimate_From_Gene_Trees/reference_approaches_source_codes/ASTRAL/ASTRAL2/ASTRAL-master-4.10.12/Astral/astral.4.10.12.jar'
                #astralcmd = 'java -jar ' + astralexec + ' -i ' + filename + ' -o ' + astraltreefile + ' > ' + outtextfile
                #os.system(astralcmd)
                
                #MakePerfectNewickTree.ConvertPerfectNewick(astraltreefile)

#---------------------------------------------------------------------------------------

if __name__=='__main__':
    main()
