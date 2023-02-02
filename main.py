import numpy as np
from Bio import SeqIO
import numpy
import matplotlib.pyplot as plt


# WARNING: Unsorted
def find_sub(seq, sub):
    # This feels like an awful way of doing this but of all the things I've tried it's the fastest by far
    start = 0
    while True:
        start = seq.find(sub, start)
        if start == -1: return
        yield start
        start += len(sub)


def site_instances(seq, site_seq):
    # Finds forwards and backwards
    instances = list(find_sub(record.seq, restriction_seq)) + list(find_sub(record.seq, restriction_seq[::-1]))
    instances.sort()
    return instances

# Read the record
record = SeqIO.read("Sequences/Long.gb", "genbank")

print("Sequence " + record.name + " of length " + str(len(record.seq)))
# Change this depending on the restriction enzyme being used
restriction_seq = "GATC"

# Count instances
instances = site_instances(record.seq, restriction_seq)

# Create empty distances array
distances = numpy.diff(instances)

print("\nInstances: " + str(len(instances)))
if len(instances) != 0:
    # Actual
    print("Mean Distance: " + str(numpy.mean(distances)))
    # Expected Values
    expected_instances = len(record.seq) / len(restriction_seq) * 1/(4 ** len(restriction_seq))
    expected_distance = len(record.seq) / expected_instances
    print("Expected Instances: " + str(expected_instances))
    print("Expected Mean Distance: " + str(expected_distance))

# Forward only log linear   Forward only log log
# Both log linear           Both log log
fig, (b_loglin, b_loglog) = plt.subplots(2)
counts, bins = np.histogram(distances, bins=100)
b_loglin.stairs(counts, bins)
b_loglog.stairs(counts, bins)
b_loglin.set_yscale("log")
b_loglog.set_yscale("log")
b_loglog.set_xscale("log")
b_loglin.set_ylabel("Frequency")
b_loglog.set_ylabel("Frequency")
b_loglin.set_xlabel("Distances")
b_loglog.set_xlabel("Distances")
# plt.hist(bins[:-1], bins, weights=counts)
plt.show()