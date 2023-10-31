                
import os
import pandas as pd
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
from kmaMat2vcf import *

read_matrix("testing1.mat", outfile="testingout.tsv")
mutdicts = getMutations(dataFile="testingout.tsv")
mutsdict2df(mutdicts, fout="testingFull.csv", sampleName="TEST")
