import numpy as np

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