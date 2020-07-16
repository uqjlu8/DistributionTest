# -*- coding: utf-8 -*-
"""
Created on Thu Jul 16 21:02:51 2020

@author: s4142554

#creates matrices of beta values from methylation files in InputFolder
"""



from operator import itemgetter, is_not
from collections import OrderedDict
import datetime
import time
import sys
from math import log10, sqrt
import re
from functools import partial 
from time import gmtime, strftime
import os
import gc
import shutil
import numpy as np
import pandas as pd
import random

from distributiontest import MatrixDistributionTest
###############################################################################
script = os.path.basename(__file__).split('.',1)[0]
script_dir = os.path.dirname(__file__) #<-- absolute dir the script is in
###############################################################################
#change to folder where all beta files are kept.
InputFolder = "C:/Users/s4142554/Desktop/mDist_014/mDist_014_Input/" 

ResultsFolder = "C:/Users/s4142554/Desktop/mDist_014/mDist_014_Results/"
###############################################################################

resultMatrix = os.path.join(script, ResultsFolder, "testMatrix.txt") #testMatrix.txt can be changed

matrixDict = {"ProbeID": []}

#get list of files in inputfolder
BetaFiles = os.listdir(InputFolder)

print ("Reading probes into dictionary...")

#go through each file individually
for i, betafile in enumerate(BetaFiles):
    betafile = os.path.join(script, InputFolder, betafile)
    
    #add name to matrix dict 
    name = "_".join(["sample", str(i+1)]) #name = sample_1, sample2 etc. #can be replaced with another naming system
    #add name to matrixDict
    matrixDict["ProbeID"].append(name)
    
    #open file and go through each probe and beta value
    with open(betafile, 'r') as file:
        for line in file:
            line = line.split()
            if line[0].startswith("Composite") == False:
                probe = line[0]
                beta = line[1]
                
                if probe in matrixDict:
                    matrixDict[probe].append(beta)
                else:
                    #check if NA needs to be added
                    if i != 0:
                        matrixDict[probe] = ["NA" for n in range(i)]
                        matrixDict[probe].append(beta)
                    elif i == 0:
                        matrixDict[probe] = [beta]
                        
    #check if NA needs to be added to fill in any gaps in the dict
    for probe in matrixDict:
        if len(matrixDict[probe]) != i+1:
            matrixDict[probe].append("NA")

print ("Writing matrix to file...")
with open(resultMatrix, 'w') as rfile:
    for probe in matrixDict:
        rline = [probe]
        rline.extend(matrixDict[probe])
        rline = "\t".join(rline)+"\n"
        rfile.write(rline)
                     
                    
        