import os
from Bio import SeqIO
from insights import *
import numpy as np
import matplotlib.pyplot as plt


# WARNING: Unsorted
def find_sub(str, sub):
    """Returns all indices of sub in str"""
    # This feels like an awful way of doing this but of all the things I've tried it's the fastest by far
    start = 0
    while True:
        start = str.find(sub, start)
        if start == -1: return
        yield start
        start += len(sub)


def site_instances(seq, site_seq):
    """Returns all indices of a restriction site sequence in seq. Only use with valid palindromic restriction site
    sequences"""
    # Concatenate forwards and backwards searches, and then sort
    instances = list(find_sub(seq, site_seq)) + list(find_sub(seq, site_seq[::-1]))
    instances.sort()
    return instances


# paths = os.listdir("Sequences/Human Genome T2T/")
paths = ["CP068271.2.gb"]

for path in paths:

    # Read the record
    record = SeqIO.read("Sequences/Human Genome T2T/" + path, "genbank")
    print("Sequence " + record.name + " of length " + str(len(record.seq)))

    # Four cutters: DpnII and MboI -- GATC
    # Six cutters: Ncd -- CCATGG, HindIII -- AAGCTT
    restriction_seq = "GATC"
    seq = record.seq

    # Remove unwanted region from sequence
    seq = seq[int(2e7):int(5.5e7)] + seq[int(6.5e7):-int(2e7)]

    instances = site_instances(seq, restriction_seq)
    distances = np.diff(instances)

    print("\nInstances: " + str(len(instances)))
    print("Mean Distance: " + str(np.mean(distances)))
    print("Cutting Probability (sites / sequence length): " + str(len(instances) / len(seq)))

    # make the figure
    fig, axs = plt.subplots(2, 2)
    fig.tight_layout(pad=3.0)

    site_distribution(axs[0][0], instances)
    # avg_distance_between_sites(axs[0][1], instances)
    shifted_return_density(axs[0][1], distances, 1)
    # autocorrelation(axs[1][0], distances)
    shifted_returns(axs[1, 0], distances, 1)
    distances_scatter(axs[1][1], distances, "linear", "log")


    def exp(x, a, b):
        return a * np.exp(-b * x)

    popt, pcov = fit_distances(axs[1][1], exp, distances)

    print("Cutting Probability (best fit line): " + str(1 - np.exp(-popt[1])))
    print("Characteristic length: " + str(1/popt[1]))
    print("Exponential Fit: " + str(round(popt[0], 3)) + " * e^(-" + str(popt[1]) + "x)")# + " + str(popt[2]))
    print("Covariance: " + str(pcov))

    fig.suptitle(record.name + " with restriction sequence " + restriction_seq)
    plt.show()
    # plt.savefig("C:/Users/lbjun/OneDrive/Documents/School/Di Pierro Lab/Fig Output/Good Sequences/Ncd/" + record.name + ".png")
    plt.close(fig)