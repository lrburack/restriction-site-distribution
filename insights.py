import numpy as np
import statsmodels.api as sm
from scipy.optimize import curve_fit


def site_distribution(ax, instances, bincount=200):
    """Creates a histogram restriction site frequency along a chromosome"""
    counts, bins = np.histogram(instances, bins=bincount)
    ax.plot(bins[1:], counts)
    ax.set_xlabel("Position (bp)")
    ax.set_ylabel("Site frequency (binned)")
    ax.set_title("Site positions")


def avg_distance_between_sites(ax, instances, bincount=200):
    """Creates a plot showing the average distance between sites at locations along a chromosome"""
    counts, bins = np.histogram(instances, bins=bincount)
    # Find bin size
    binsize = bins[1] - bins[0]
    # Remove bins with zero instances
    bins = bins[1:][counts != 0]
    counts = counts[counts != 0]
    ax.plot(bins, binsize / counts)
    ax.set_xlabel("Position (bp)")
    ax.set_ylabel("Average site distance (binned)")
    ax.set_title("Site distances")


def distances_scatter(ax, distances, xscale="linear", yscale="linear", bincount=300):
    """Creates a histogram restriction site frequency along a chromosome"""
    counts, bins = np.histogram(distances, bins=bincount)
    x = bins[:-1]  # Use the low end of each bin as the x value

    ax.scatter(x, counts, s=5, color="lightcoral")
    ax.set_xscale(xscale)
    ax.set_yscale(yscale)
    ax.set_xlabel("Distance between sites (bp)")
    ax.set_ylabel("Frequency (binned)")
    ax.set_title("Distance distribution")


def autocorrelation(ax, distances, nlags=100):
    acorr = sm.tsa.acf(distances, nlags=10000)
    ax.plot(range(len(acorr)), acorr)
    ax.set_title("ACF")
    ax.set_yscale("log")
    ax.set_xscale("log")
    ax.set_ylabel("ACF")
    ax.set_xlabel("Lag")


def fit_distances(ax, func, distances, bins=200):
    counts, bins = np.histogram(distances, bins=bins)
    x = bins[:-1]
    popt, pcov = curve_fit(func, x, counts, maxfev=60000)
    y = func(x, *popt)
    ax.plot(x, y, label=func.__name__ + " fit", color="navy", linestyle="dashed")
    return popt, pcov


def shifted_returns(ax, distances, lag):
    if lag == 0:
        ax.scatter(distances, distances, s=5, color="lightcoral")
    else:
        ax.scatter(distances[lag:], distances[:-lag], s=5, color="lightcoral")
    # ax.scatter(distances[lag:], distances[:-lag], s=5, color="lightcoral")
    ax.set_ylabel("Shifted up " + str(lag))
    ax.set_xlabel("Shifted back " + str(lag))


def shifted_return_density(ax, distances, lag=1, bincount=300):
    # Get the distance from the origin of each point in the shift. Bin these distances and plot them.
    if lag == 0:
        shifted_distances = np.sqrt(np.square(distances) + np.square(distances))
    else:
        shifted_distances = np.sqrt(np.square(distances[lag:]) + np.square(distances[:-lag]))

    counts, bins = np.histogram(shifted_distances, bincount)
    ax.plot(bins[1:], counts)
    ax.set_ylabel("Density")
    ax.set_xlabel("Shift dist from origin")


# def regions_of_interest(instances, distances, z_thresh=10, consecutive_thresh=10):
#     """Identify ranges in the sequence containing elevated or depressed site frequencies"""
#     mean = np.mean(distances)
#     stdev = np.std(distances)
#
#     problematic_inds = np.arange(0, len(distances))[distances > mean + z_thresh * stdev]
#
#     ranges = np.split(problematic_inds, np.where(np.diff(problematic_inds) != 1)[0] + 1)
#
#     return problematic_inds