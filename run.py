# -*- coding: utf-8 -*-
"""
Created on Mon Jun 29 14:09:14 2020

@author: s4142554
"""
from collections import OrderedDict
import getopt
import sys
import os
import argparse
from distributiontest import MatrixDistributionTest
###############################################################################
script = os.path.basename(__file__).split('.',1)[0]
script_dir = os.path.dirname(__file__) #<-- absolute dir the script is in
###############################################################################
#replace InputMatrixFile with matrix of your choice
InputMatrixFile = 'TCGA-BRCA_TissueNormal_merged-matrix.txt'
ResultFIle = 'distribution_results_summary.txt'

if os.path.exists(InputMatrixFile) == True:
    try:
        MatrixDistributionTest(InputMatrixFile, ResultFIle)
    except  MemoryError as err:
        print (err)
        