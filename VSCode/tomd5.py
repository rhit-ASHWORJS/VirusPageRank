import csv
import os
# import numpy as np
# import matplotlib.pyplot as plt
# import matplotlib as mpl
# from matplotlib import cm
# # from colorspacious import cspace_converter
# from collections import OrderedDict
from tkinter.filedialog import askopenfilename

# import ExperimentStorer
# import GraphMaker
# import statistics

#csvsize workaround---------------------------------------
import sys
import csv
maxInt = sys.maxsize

while True:
    # decrease the maxInt value by factor 10 
    # as long as the OverflowError occurs.

    try:
        csv.field_size_limit(maxInt)
        break
    except OverflowError:
        maxInt = int(maxInt/10)
#-------------------------------------------------------

filename = askopenfilename()
f = open("output.txt", "w")
numtowrite = 100
numwritten = 0

with open(filename) as csvfile:
    datareader = csv.reader(csvfile, delimiter='\t')
    for row in datareader:
        f.write(row[0] + "\n")
        numwritten = numwritten + 1
        if(numwritten == numtowrite):
            break

f.close()
