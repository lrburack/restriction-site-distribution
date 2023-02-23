import math
from Bio import SeqIO
from statsmodels.graphics.tsaplots import plot_acf

from insights import *
from processSequences import *
import numpy as np
import matplotlib.pyplot as plt

# Supress warnings
import warnings
warnings.filterwarnings("ignore")


def remove_outside(str, ranges):
    result = ''
    for range in ranges:
        result += str[range[0]:range[1]]
    return result


# This seems silly but it is actually saving time I promise
paths = ["CM000663.1.gb", "CM000669.1.gb"]
names = ["Chr1", "Chr7"]
# keep_region = [[int(1.2e8), int(1.6e8)], [int(.5e8), int(.7e8)]]
keep_region = [
    [[int(.2e8), int(1.2e8)], [int(1.5e8), -int(.2e8)]],
    [[int(.2e8), int(.5e8)], [int(.7e8), -int(.2e8)]]
]

restriction_enzymes = ["MboI", "NcoI", "HindIII"]
restriction_sequences = ["GATC", "CCATGG", "AAGCTT"]


# For each chromosome enzyme pair we want to create some figure
for path_ind in range(len(paths)):
    record = SeqIO.read("Sequences/Hg19 Assembly for Paper/" + paths[path_ind], "genbank")
    seq = remove_outside(record.seq, keep_region[path_ind])
    print("\n\nFilename: " + names[path_ind])
    print("Keeping ranges (bp): " + str(keep_region[path_ind]))

    for restr_ind in range(len(restriction_sequences)):
        instances = site_instances(seq, restriction_sequences[restr_ind])
        distances = np.diff(instances)

        # make the figure
        fig, ax = plt.subplots(1)

        fig.tight_layout(pad=2)

        # The right bin size is about 100 * average length
        # bincount = len(record.seq) // 100000
        # # Round to the nearest power of ten
        # # bincount = np.power(10, math.floor(np.log10(bincount)))
        # val, unit = bp_notation(len(seq) // bincount)
        #
        # # site distribution figure
        # site_distribution(ax, instances, bincount=bincount)
        # ax.set_ylabel(restriction_enzymes[restr_ind] + " sites per " + '%g'%val + unit)

        # # difference in distance per bin = 100bp
        # bincount = int(max(distances) // 100)
        # distances_scatter(ax, distances, "linear", "log", bincount=bincount)
        #
        # def exp(x, a, b):
        #     return a * np.exp(-b * x)
        # popt, pcov, r_squared = fit_distances(ax, exp, distances, bincount)
        #
        # # Figure title
        # fig.suptitle(names[path_ind] + " with " + restriction_enzymes[restr_ind])
        #
        # plt.savefig(
        #     "C:/Users/lbjun/OneDrive/Documents/School/Di Pierro Lab/Restriction Site Distribution Research/Fig Output/Sliced "
        #     + names[path_ind] + " " + restriction_enzymes[restr_ind])

        # poisson_region_length = 5 * 4 ** len(restriction_sequences[restr_ind])
        # poisson_region_count = len(seq) // poisson_region_length
        #
        # instances = site_instances(seq, restriction_sequences[restr_ind])
        # distances = np.diff(instances)
        #
        # popt, pcov, r_squared = poisson_distribution(ax, instances, regioncount=poisson_region_count)
        # print(r_squared)
        print("plotting")
        plot_acf(distances, ax=ax, lags=range(0, 1000, 10))

        print("Restriction enzyme: " + restriction_enzymes[restr_ind])
        print("- Restriction site: " + restriction_sequences[restr_ind])
        print("- Number of instances: " + str(len(instances)))
        print("- Average distance (bp): " + str(len(seq) // len(instances)))
        # print("- Exponential Fit: " + str(popt[0]) + " * e^(-" + str(popt[1]) + "* x)")
        # print("Optimized Params: " + str(popt))
        # print("- r^2: " + str(r_squared))

        plt.savefig(
            "C:/Users/lbjun/OneDrive/Documents/School/Di Pierro Lab/Restriction Site Distribution Research/Fig Output/"
            + names[path_ind] + " " + restriction_enzymes[restr_ind])
