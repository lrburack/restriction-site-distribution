import math

import numpy as np
from Bio import SeqIO
import numpy
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit


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

def view_sites(seq, site_seq):
    instances = site_instances(seq, site_seq)
    distances = numpy.diff(instances)

    print("\nInstances: " + str(len(instances)))
    if len(instances) == 0:
        return
    # Actual
    print("Mean Distance: " + str(numpy.mean(distances)))
    # Expected Values
    expected_instances = len(seq) / len(restriction_seq) * 1 / (4 ** len(restriction_seq))
    expected_distance = len(seq) / expected_instances
    print("Expected Instances: " + str(expected_instances))
    print("Expected Mean Distance: " + str(expected_distance))

    # Remove large outliers
    remove_above_zscore = 10
    stdev = numpy.std(distances)
    mean = numpy.mean(distances)
    print("\nRemoved " + str(sum(distances > mean + stdev * remove_above_zscore)) + " outliers.")
    desert_region_percentage = round(
        100 * numpy.sum(distances[distances > mean + stdev * remove_above_zscore]) / len(seq), 5)
    print("These outliers (desert regions) compose " + str(desert_region_percentage) + "% of the total sequence")
    desert_site_inds = numpy.where(distances > mean + stdev * remove_above_zscore)[0]
    print("Desert inds: " + str(desert_site_inds))
    print("Start: " + str(instances[desert_site_inds[0]]) + "\tEnd: " + str(instances[desert_site_inds[-1]]) + "\tLength: " + str(instances[desert_site_inds[-1]] - instances[desert_site_inds[0]]))
    distances = distances[distances < mean + stdev * remove_above_zscore]

    # make the figure
    fig, ((linlin, linlog), (loglog, exp_err)) = plt.subplots(2, 2)
    fig.tight_layout(pad=3.0)

    # Bin counts
    counts, bins = np.histogram(distances, bins=200)
    # Dont include zeros because they cause problems with logs and also look bad
    mask = counts != 0
    print("Removed " + str(sum(mask == 0)) + " empty bins.")
    x = bins[:-1][mask]
    y = counts[mask]

    linlin.set_title("Lin Lin")
    linlog.set_title("Lin Log")
    loglog.set_title("Log Log")
    exp_err.set_title("Exponential Fit Error")

    # Calculate exponential best fit line
    exp_a, exp_b = numpy.polyfit(x, numpy.log(y), 1)
    # Calculate power law best fit
    pow_a, pow_b = numpy.polyfit(numpy.log(x), numpy.log(y), 1)

    # Plot data
    axs = [linlin, linlog, loglog]
    for ax in axs:
        ax.scatter(x, y, s=5, color="lightcoral")
        ax.plot(x, numpy.exp(exp_b) * numpy.exp(exp_a * x), color="navy", linestyle="dashed", label="Exponential Fit")
        ax.plot(x, numpy.power(x, pow_a) * numpy.exp(pow_b), color="royalblue", linestyle="dotted", label="Power Fit")
        ax.set_ylabel("Frequency")
        ax.set_xlabel("Distances (bp)")

    ymax = numpy.max(y)
    # Additional formatting stuff
    linlin.legend()
    loglog.set_yscale("log")
    loglog.set_xscale("log")
    linlog.set_yscale("log")
    linlog.set_ylim([1, ymax])
    loglog.set_ylim([1, ymax])
    linlin.set_ylim([0, ymax])

    err = y - numpy.exp(exp_b) * numpy.exp(exp_a * x)
    stdev = numpy.std(err)
    exp_err.scatter(x, err/stdev, s=5, color="lightcoral")
    exp_err.plot([x[0], x[-1]], [0, 0], color="black", linestyle="dashed", label="Zero")
    exp_err.set_ylabel("Z Score")
    exp_err.set_xlabel("Distances (bp)")
    exp_err.legend()

    return fig

# Read the record
record = SeqIO.read("Sequences/CM039034.1.gb", "genbank")
print("Sequence " + record.name + " of length " + str(len(record.seq)))

# Four cutters: DpnII and Mbd -- GATC
# Six cutters: Ncd -- CCATGG, HindIII -- AAGCTT
restriction_seq = "CCATGG"

fig = view_sites(record.seq, restriction_seq)
fig.suptitle(record.name + " with restriction sequence " + restriction_seq)
plt.show()