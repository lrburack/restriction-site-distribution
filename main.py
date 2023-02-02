import numpy as np
from Bio import SeqIO
import numpy
import matplotlib.pyplot as plt


def find_all(a_str, sub):
    start = 0
    while True:
        start = a_str.find(sub, start)
        if start == -1: return
        yield start
        start += len(sub)

# Read the record
record = SeqIO.read("Sequences/Long.gb", "genbank")
print("Sequence " + record.name + " of length " + str(len(record.seq)))
# Change this depending on the restriction enzyme being used
restriction_seq = "GATC"

# Count instances
instances = list(find_all(record.seq, restriction_seq))

# Create empty distances array
distances = [None] * (len(instances) - 1)
# Find distances
for i in range(1, len(instances)):
    distances[i-1] = instances[i] - instances[i-1]

print("Instances: " + str(len(instances)))
if len(instances) != 0:
    # Actual
    # print("Distances: " + str(distances))
    print("Mean Distance: " + str(numpy.mean(distances)))

    # Expected Values
    expected_instances = len(record.seq) / len(restriction_seq) * 1/(4 ** len(restriction_seq))
    expected_distance = len(record.seq) / expected_instances
    print("Expected Instances: " + str(expected_instances))
    print("Expected Mean Distance: " + str(expected_distance))

counts, bins = np.histogram(distances, bins=100)
plt.stairs(counts, bins)
print("Bins: " + str(len(bins)))
print("Counts: " + str(len(counts)))
plt.yscale("log")
plt.ylabel("Frequency")
plt.xlabel("Distances")
# plt.hist(bins[:-1], bins, weights=counts)
plt.show()