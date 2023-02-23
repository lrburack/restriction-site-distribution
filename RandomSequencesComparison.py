import random
from insights import *
from processSequences import *
import matplotlib.pyplot as plt
from statsmodels.graphics.tsaplots import plot_acf
from Bio import SeqIO
from timeit import timeit


def randomSeq(length=int(), letters="GCTA"):
    return''.join(random.choices(letters, k=length))


restriction_seq = "CCATGG"
sequence_length = int(1e7)

instances = site_instances(randomSeq(sequence_length), restriction_seq)
distances = np.diff(instances)

fig, ax = plt.subplots(1)
fig.tight_layout(pad=2)

plot_acf(distances, ax=ax, lags=range(0, 1000, 10))

plt.show()