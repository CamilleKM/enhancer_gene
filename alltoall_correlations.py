#Calcul the Pearson correlation coefficient between the promoter and
#all DHS peak within less than 500kb.
#Usage:
#python alltoall_correlations.py [-p prom_file.bed] [-d 500000] DNase_DHS.bed
#CPU = 16 by default

import argparse as ap
from math import log10
import numpy as np
import sys
from multiprocessing import Pool, sharedctypes
import warnings


class Dhs(object): #For each DHS.
    """
    Define each DHS peak as an object and compute firsts calculs: centered
    logcounts and standard deviation.
    """
    def __init__(self,chrom,start,end,name,counts):
        self._chrom = chrom
        self._start = int(start)
        self._end = int(end)
        self._name = name
        self.counts = map(float,counts) #All counts are floats.
        self.sigma = np.std(self.log_counts)
        self.moy = np.mean(self.log_counts)
        self.log_counts_prime = np.array([i-self.moy for i in self.log_counts])

    def __str__(self):
        return ("%s\t%s\t%s\t%s" % (self.chrom,self.start,self.end,self.name))

    def __repr__(self):
        return "%s:%s-%s\t%s\t%s" % (self._chrom, self._start, self._end,
                                     self._id, floats_to_str(self.counts))

    @property
    def chrom(self):
        return self._chrom

    @property
    def start(self):
        return self._start

    @property
    def end(self):
        return self._end

    @property
    def name(self):
        return self._name

    @property
    def length(self):
        """
        Return the length of the counts vector.
        """
        return (len(self.counts))

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

    def within_neighbourhood(self,dhs_interval,distance):
        """
        Return True if DHS interval is within 500kb of self (the promoter)
        Param dhs_interval: vector [chrom,start,end] and distance (interger)
        Return: bolean value.
        """
        if self.chrom == dhs_interval[0] and (self.end < dhs_interval[1] and self.end+distance >= dhs_interval[1]):
            return True
        else:
            return False


def read_dhs_file(dhs_file):
    """
    Simply read the DHS file into an array of Dhs objectsfields
    """
    all_dhs = []  # Vector with all the object DHS from file_dnase.
    with open(dhs_file, "r") as f:
        for line in f:
            fields = line.rstrip().split("\t")
            if fields[0] != "#chr":
                chrom, start, end, identifier = fields[0:4]
                start = int(start)
                end = int(end)
                counts = [float(c) for c in fields[4:]]
                all_dhs.append(Dhs(chrom, start, end, identifier, counts))
    return all_dhs


def write_correlations(output):
    header = "\t".join(["#chr", "start", "end", "ID",
                        "chr", "start", "end", "ID",
                        "Pearson_corr"])
    print(header)
    for dhs1,dhs2,corr in output:
        print("%s\t%s\t%5.4f" % (dhs1, dhs2, corr))


def dhs_arraytonumpy (all_dhs):
    """
    Convert values from the class Dhs (logcounts and std) to numpy arrays for
    multiprocessing.
    Param: the table with all the object.
    Return: two numpy array.
    """
    table_logcounts = []
    table_std = []
    for i in all_dhs:
        table_logcounts.append(i.log_counts_prime)
        table_std.append(i.sigma)
    return (np.array(table_logcounts), np.array(table_std))

def pearson_correlation (logcount1,logcount2,sigma1,sigma2):
    """
    Compute the pearson correlation coefficient by hand thanks to calculs made
    in the class Dhs.
    Params: take two Dhs objects containing all the precomputed values.
    Return: pearson correlation coefficient.
    """
    numerator = np.dot(logcount1,logcount2)
    coefP = numerator/(10*sigma1*sigma2)
    return coefP


def parallel_correlations(args):
    """
    Compute all the correlation coefficient calcul of the Dhs in start from
    start+1 to end.
    Param: to_compute with the windows border index for each Dhs.
    Return: list of the correlation coefficient.
    """
    start,end = args

    with warnings.catch_warnings(): #Warning due to the use of global variables.
        warnings.simplefilter('ignore', RuntimeWarning)
        logcounts = np.ctypeslib.as_array(logcounts_global)
        sigma = np.ctypeslib.as_array(sigma_global)

    correlations = []
    for j in range(start+1, end+1):
        r = pearson_correlation(logcounts[start],logcounts[j],sigma[start],sigma[j])
        correlations.append(((start, j), r))
    return correlations


def _init_parallel(logcounts_to_populate,sigma_to_populate):
    """
    Each pool process calls this initializer.
    Load the arrays to be populated into that process's global namespace
    """
    global logcounts_global
    global sigma_global
    logcounts_global = logcounts_to_populate
    sigma_global = sigma_to_populate


def parallel_implementation_corr(logcounts, sigma, to_compute, num_proc=1):
    """
    Create shared variables and pool for multiprocessing.
    Params: logcounts and sigma numpy array, to_compute list with the index and
    the number of processeurs to use.
    Return: all the correlations and index
    """
    # Define the corresponding shared ctype arrays
    tmp_logcounts = np.ctypeslib.as_ctypes(logcounts) #Create ctypes object from numpy array
    shared_logcounts = sharedctypes.Array(tmp_logcounts._type_, tmp_logcounts, lock=False)
    #                       Return a ctypes array allocated from shared memory.
    #                       tmp_dhs._type_: type of the returned array elements
    #                       tmp_dhs: initialize the array
    #                       lock=True (default): synchronize access to the value
    #                       lock=False: access to returned object not protected

    tmp_sigma = np.ctypeslib.as_ctypes(sigma)
    shared_sigma = sharedctypes.Array(tmp_sigma._type_,tmp_sigma, lock=False)

    pool = Pool(processes=num_proc, initializer=_init_parallel,
                initargs=(shared_logcounts, shared_sigma, ))
    #                       controls a pool of worker processes
    #                       processes: number of worker processes to use
    #                       initializer is not None: each worker process will
    #                       call initializer(*initargs)
    result = pool.map(parallel_correlations, to_compute)
    #                       supports only one iterable argument,
    #                       blocks until the result is ready

    results = [item for sublist in result for item in sublist]
    return results


def sort_correlations(correlations,all_dhs):
    """
    Sort the result from multiprocessing and add the data corresponding to each
    coefficient.
    Params: list with index and correlation coefficient and list of all the
    object Dhs
    Result: sorted numpy array with the results to print.
    """
    sorted_cor = [[all_dhs[cell[0][0]],all_dhs[cell[0][1]],cell[1]] \
                  for cell in sorted(correlations,key=lambda x:(x[0][0], x[0][1]))]
    return np.array(sorted_cor)


if __name__ == '__main__':
    #Input from command line.
    parser = ap.ArgumentParser()
    parser.add_argument("file_dnase", help="file with the read count per peak")
    parser.add_argument("-p", "--with_promoter", help="if active, give a file \
                        with the read count per peak for each promoter")
    parser.add_argument("-d", "--distance", help="distance maximum between peaks\
                        to calcul the correlation. (default:500000)",
                        default=500000,type=int)
    args = parser.parse_args()

    all_dhs = read_dhs_file(args.file_dnase)

    to_compute = []
    if args.with_promoter:
        print("Avec le promoteur c'est plus dur, option a venir")
    else:
        for i in range(len(all_dhs)-1):
            right_border = i+1
            while all_dhs[i].within_neighbourhood(\
                                                  all_dhs[right_border].interval,\
                                                  args.distance):
                right_border += 1
                if right_border == len(all_dhs):
                    break
            to_compute.append([i,right_border-1]) #For each peak, we compute
            #correlation with all the peaks in this window.

    num_proc = 16
    (logcounts,sigma) = dhs_arraytonumpy(all_dhs)
    output = parallel_implementation_corr(logcounts, sigma, to_compute, num_proc)
    sorted_correlations = sort_correlations(output,all_dhs)

    write_correlations(sorted_correlations)
