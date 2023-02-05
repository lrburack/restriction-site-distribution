import numpy as np
from scipy.optimize import curve_fit

def site_distribution(ax, instances, bins=200):
    """Creates a histogram restriction site frequency along a chromosome"""
    counts, bins = np.histogram(instances, bins=bins)
    ax.stairs(counts, bins)
    ax.set_xlabel("Position (bp)")
    ax.set_ylabel("Site frequency (binned)")
    ax.set_title("Site positions")


def avg_distance_between_sites(ax, instances, bins=200):
    """Creates a scatter showing the average distance between sites at locations along a chromosome"""
    counts, bins = np.histogram(instances, bins=bins)
    # Find bin size
    binsize = bins[1] - bins[0]
    # Remove bins with zero instances
    bins = bins[1:][counts != 0]
    counts = counts[counts != 0]
    ax.plot(bins, binsize / counts)
    ax.set_xlabel("Position (bp)")
    ax.set_ylabel("Average site distance (binned)")
    ax.set_title("Site distances")


def distances_scatter(ax, distances, xscale="linear", yscale="linear", bins=300):
    """Creates a histogram restriction site frequency along a chromosome"""
    counts, bins = np.histogram(distances, bins=bins)
    x = bins[:-1]  # Use the low end of each bin as the x value
    ax.scatter(x, counts, s=5, color="lightcoral")
    ax.set_xscale(xscale)
    ax.set_yscale(yscale)
    ax.set_xlabel("Distance between sites (bp)")
    ax.set_ylabel("Frequency (binned)")
    ax.set_title("Distance distribution")




# def regions_of_interest(instaces, distances, z_thresh=10):
#     mean = np.mean(distances)