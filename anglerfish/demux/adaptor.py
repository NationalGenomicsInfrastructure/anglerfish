import importlib
import os

import yaml


class Adaptor:
    def __init__(self, adaptors, delim, adaptor, i7_index=None, i5_index=None):
        self.i7 = AdaptorPart(adaptors[adaptor]["i7"], adaptor, delim, i7_index)
        self.i5 = AdaptorPart(adaptors[adaptor]["i5"], adaptor, delim, i5_index)
        self.name = f"{adaptor}"
        self.delim = delim

    def get_i7_mask(self, insert_Ns=True):
        return self.i7.get_mask(insert_Ns)

    def get_i5_mask(self, insert_Ns=True):
        return self.i5.get_mask(insert_Ns)

    def get_fastastring(self, insert_Ns=True):
        fasta_i5 = f">{self.name}_i5\n{self.get_i5_mask(insert_Ns)}\n"
        fasta_i7 = f">{self.name}_i7\n{self.get_i7_mask(insert_Ns)}\n"
        return fasta_i5 + fasta_i7


class AdaptorPart:
    def __init__(self, sequence, name, delim, index):
        self.sequence = sequence
        self.name = name
        self.delim = delim
        self.index = index

    def has_index(self):
        return self.sequence.find(self.delim) > -1

    def len_before_index(self):
        return self.sequence.find(self.delim)

    def len_after_index(self):
        return len(self.sequence) - self.sequence.find(self.delim) - len(self.delim)

    def get_mask(self, insert_Ns):
        if self.has_index():
            if not insert_Ns:
                return self.sequence.replace(self.delim, "")
            else:
                return self.sequence.replace(self.delim, "N" * len(self.index))
        else:
            return self.sequence


# Fetch all adaptors
def load_adaptors(raw=False):
    p = importlib.resources.files("anglerfish.config").joinpath("adaptors.yaml")
    assert isinstance(p, os.PathLike)

    adaptors_raw = []
    with open(p) as f:
        adaptors_raw = yaml.safe_load(f)

    if raw:
        return adaptors_raw
    adaptors = []
    for adaptor in adaptors_raw:
        adaptors.append(
            Adaptor(adaptors_raw, "-NNN-", adaptor, i7_index=None, i5_index=None)
        )

    return adaptors
