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


paths = os.listdir("Sequences/Human Genome")
for path in paths:

    # Read the record
    record = SeqIO.read("Sequences/Human Genome/" + path, "genbank")
    print("Sequence " + record.name + " of length " + str(len(record.seq)))

    # Four cutters: DpnII and MboI -- GATC
    # Six cutters: Ncd -- CCATGG, HindIII -- AAGCTT
    restriction_seq = "GATC"
    seq = record.seq

    instances = site_instances(seq, restriction_seq)
    distances = np.diff(instances)

    print("\nInstances: " + str(len(instances)))
    print("Mean Distance: " + str(np.mean(distances)))

    # Remove large outliers
    # remove_above_zscore = 10
    # stdev = np.std(distances)
    # mean = np.mean(distances)
    # print("\nRemoved " + str(sum(distances > mean + stdev * remove_above_zscore)) + " outliers.")
    # print("Distances that were removed: " + str(distances[distances > mean + stdev * remove_above_zscore]))
    # # Preprocess
    # distances = distances[distances < mean + stdev * remove_above_zscore]


    # make the figure
    fig, axs = plt.subplots(2,1)
    fig.tight_layout(pad=3.0)

    site_distribution(axs[0], instances)
    avg_distance_between_sites(axs[1], instances)

    fig.suptitle(record.name + " with restriction sequence " + restriction_seq)
    plt.savefig("C:/Users/lbjun/OneDrive/Documents/School/Di Pierro Lab/Fig Output/" + record.name + ".png")
    plt.close(fig)