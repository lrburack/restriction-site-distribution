from Bio import SeqIO
import numpy as np
import matplotlib.pyplot as plt


# This is unreadable and terrible, but I just want it done
def base_histogram(seq, char, bincount):
    binlength = len(seq) // bincount
    bins = np.arange(0, len(seq), binlength)
    counts = [None] * bincount

    for i in range(1, len(bins)):
        counts[i - 1] = seq[bins[i - 1]:bins[i]].count(char)

    return bins, counts, binlength


bincount = 10000
path = "CP068271.2.gb"
record = SeqIO.read("Sequences/Human Genome T2T/" + path, "genbank")
seq = str(record.seq)

fig, [base_frequency, amine_freq] = plt.subplots(2)

bins, a_counts, binlength = base_histogram(seq, "A", bincount)
at_counts = np.add(a_counts, base_histogram(seq, "T", bincount)[1])
at_normcounts = np.divide(at_counts, binlength)
base_frequency.plot(bins[:-1], at_normcounts, label="AT")

gc_counts = np.add(base_histogram(seq, "G", bincount)[1], base_histogram(seq, "C", bincount)[1])
gc_normcounts = np.divide(gc_counts, binlength)
base_frequency.plot(bins[:-1], gc_normcounts, label="GC")

together_counts = np.add(at_normcounts, gc_normcounts)
base_frequency.plot(bins[:-1], together_counts, label="Both")

base_frequency.plot([bins[0], bins[-1]], [.5, .5], color="black", linestyle="dashed")
base_frequency.legend()
base_frequency.set_ylim([0, 1.1])
base_frequency.set_title("Base Frequency")
base_frequency.set_ylabel("Frequency")
base_frequency.set_xlabel("Position (bp)")

amine_freq.plot(bins[:-1], np.add(np.multiply(at_counts, .5), gc_counts))

plt.show()


