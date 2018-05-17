#Calcul the Pearson correlation coefficient between the promoter and all DHS peak
#within less than 500kb.
#Usage: python Prom-other_correlation_eco.py DNase_DHS.bed DHS_promoter.bed > prom-other_correlation.txt

import argparse as ap
import scipy.stats
from math import log10
import numpy as np
import sys

class DHS(object): #For each DHS.
    """
    Define each DHS peak as an object and compute firsts calculs: centered
    logcounts and standard deviation.
    """
    def __init__(self,fields):
        self.chrom = fields[0]
        self.start = int(fields[1])
        self.end = int(fields[2])
        self.ID = fields[3]
        self.counts = map(float,fields[4:]) #All counts are floats.
        self.sigma = np.std(self.log_counts)
        self.moy = np.mean(self.log_counts)
        self.log_counts_prime = np.array([i-self.moy for i in self.log_counts])

    def __str__(self):
        return ("%s\t%s\t%s\t%s" % (self.chrom,self.start,self.end,self.ID))

    @property
    def log_counts(self):
        """
        Avoid false correlation and NA values
        Return: vector with the log counts.
        """
        epsilon = 0.001
        return [log10(i+epsilon) for i in self.counts]

    @property
    def interval(self):
        """
        Create a vector about the position of the DHS peak.
        Return: vector [chrom,start,end]
        """
        return [self.chrom,self.start,self.end]

    def is_within_neighbourhood(self,DHS_interval):
        """
        Return True if DHS interval is within 500kb of self (the promoter)
        Param DHS_interval: vector [chrom,start,end]
        Return: bolean value.
        """
        if ((self.chrom == DHS_interval[0])
                 and ((self.start > DHS_interval[2] and self.start <= DHS_interval[2]+500000)
                 or (self.end < DHS_interval[1] and self.end+500000 >= DHS_interval[1]))) or self.interval == DHS_interval:
            return True
        else:
            return False


def DHS_correlation (DHS1, DHS2):
    """
    Compute the pearson correlation coefficient by hand thanks to calculs made
    in the class DHS.
    Params: take two DHS objects containing all the precomputed values.
    Return: pearson correlation coefficient.
    """
    numerator = np.dot(DHS1.log_counts_prime,DHS2.log_counts_prime)
    coefP = numerator/(10*DHS1.sigma*DHS2.sigma)
    return str(coefP)


#Input from command line.
parser = ap.ArgumentParser()
parser.add_argument("file_dnase")
parser.add_argument("file_prom")
args = parser.parse_args()

all_DHS = [] #Vector with all the object DHS from file_dnase.
with open (args.file_dnase,"r") as f:
    for line in f:
        fields = line.rstrip().split("\t")
        if fields[0] != "#chr":
            all_DHS.append(DHS(fields))

left_border = 0 #Beginning of the window around the studied promoter.
output=[] #Results to print
output.append("\t".join(["#prom_chr","prom_start","prom_end","prom_ID","chr","start",
"end","ID","Pearson_corr"])) #Header line
with open (args.file_prom,"r") as f:
    dejavu=[] #Avoid to calcul the same correlation twice.
    for line in f:
        fields = line.split("\t")
        if fields[0] != "#chr":
            prom_DHS = DHS(fields)

            #while the left border is not in the neighbourhood of prom_DHS, the
            #window move forward.
            while not prom_DHS.is_within_neighbourhood(all_DHS[left_border].interval):
                left_border += 1
                if left_border > len(all_DHS):
                    break

            #Once the left border of the window is within the neighbourhood, we
            #can compute the correlation coefficient.
            i = left_border
            while prom_DHS.is_within_neighbourhood(all_DHS[i].interval):
                if left_border> len(all_DHS):
                    break
                if prom_DHS.interval != all_DHS[i].interval and [all_DHS[i],prom_DHS] not in dejavu:
                    coefP = DHS_correlation(prom_DHS,all_DHS[i])
                    output.append("\t".join([str(prom_DHS),str(all_DHS[i]),coefP,"\n"]))
                    dejavu.append([prom_DHS,all_DHS[i]])
                i += 1
                if i >= len(all_DHS):
                    break

#Print the results.
for line in output:
    print(line)
