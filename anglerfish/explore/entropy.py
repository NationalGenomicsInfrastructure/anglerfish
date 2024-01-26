import itertools
from collections import deque

import numpy as np


# Rolling window implementation
def _window(seq, n=3):
    win = deque(maxlen=n)
    for char in seq:
        win.append(char)
        if len(win) >= n:
            yield tuple(win)


def _entropy(pk, qk):
    """
    Kullback-Leibler divergence (relative entropy)
    Reference https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.entropy.html
    """
    # To avoid zeroes in pk, add a tiny fraction of observed values
    pk = pk + 1e-2

    if np.sum(pk) - 1 > 1e-5:
        pk = pk / np.sum(pk)
    return np.sum(pk * np.log(pk / qk))


# Each sequence gives matrix of position and k-mer occurence.
def _seq_to_count_matrix(dna_seq, kmer_positions, max=50, k=2):
    """max is the length of the prefix that will be considered"""
    if len(dna_seq) < max:
        max = len(dna_seq)

    matrix = np.zeros((4**k, max - k + 1))
    for i, word in enumerate(_window(dna_seq[:max], n=k)):
        matrix[kmer_positions[word], i] += 1
    return matrix


def _count_for_seqs(seqs, kmer_positions, insert_length, k=2):
    seqs = [seq for seq in seqs if len(seq) == insert_length]
    matrices = [_seq_to_count_matrix(seq, kmer_positions, max=50, k=k) for seq in seqs]
    stacked = np.stack(matrices)
    sum_a = np.sum(stacked, axis=0)
    even_dist = np.ones(4**k) / (4**k)
    entropies = [_entropy(a, even_dist) for a in sum_a.T]
    return entropies


def _extract_inserts_from_df(df):
    mim_re_cs = r"^cs:Z::[1-9][0-9]*\+([a,c,t,g]*):[1-9][0-9]*$"
    return df.cs.str.extract(mim_re_cs, expand=False).str.upper().tolist()


def calculate_relative_entropy(df_good_hits, kmer_length, insert_length):
    all_kmers = itertools.product("ACTG", repeat=kmer_length)
    kmer_positions = dict(
        (kmer, position) for kmer, position in zip(all_kmers, range(4**kmer_length))
    )
    seqs = _extract_inserts_from_df(df_good_hits)
    entropies = _count_for_seqs(seqs, kmer_positions, insert_length, k=kmer_length)
    return entropies
