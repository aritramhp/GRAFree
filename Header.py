# To parse the file from shell script
from optparse import OptionParser


# To read/ write csv file
# import itertools
# from itertools import chain
# import csv


# to call different numerical, statistical functions
# from scipy import stats
# from scipy.stats import norm
import math
import numpy
from numpy import linalg


# To run shell commands from python
import sys, os, warnings
# import shutil


# to retrieve system time
import time


# to plot data
import matplotlib
import matplotlib.pyplot as plt
# from mpl_toolkits.mplot3d import Axes3D
# from matplotlib.collections import LineCollection
# from matplotlib.colors import ListedColormap, BoundaryNorm
# import matplotlib.patches as mpatches
# import matplotlib.mlab as mlab
import gc


# Biopython package need to work with the sequence files
from Bio import SeqIO
from Bio.Alphabet import generic_dna
from Bio import Phylo
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# To Construct tree
import Bio
from Bio import SeqIO
from Bio import AlignIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.Phylo import *

import dendropy
from dendropy import TreeList, Tree, Taxon, TaxonSet, Node

# This part need to call UPGMA
from Bio import Phylo
# import Bio.Phylo.TreeConstruction
from Bio.Phylo.TreeConstruction import _DistanceMatrix
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor

# Random number generate
import random

# # Parallel processing
# import threading
# import multiprocessing
# from multiprocessing import Pool, Process, Manager
#
# # Garbage collector
# import gc
#
# # Parse URL
# from bs4 import BeautifulSoup
# from urlparse import urlparse
# import requests as req

# Save and load dictionary
import cPickle

# For deepcopy
import copy

# # For majority voting
# from collections import Counter
