import importlib
import os
import re

import yaml

idelim = re.compile(r"\<N\>")
udelim = re.compile(r"(\<U\d+\>)")
ulen = re.compile(r"\<U(\d+)\>")


class Adaptor:
    def __init__(self, adaptors, delim, adaptor_type, i7_index=None, i5_index=None):
        self.i5 = AdaptorPart(
            adaptors[adaptor_type]["i5"], adaptor_type, delim, i5_index
        )
        self.i7 = AdaptorPart(
            adaptors[adaptor_type]["i7"], adaptor_type, delim, i7_index
        )
        self.i5_index = i5_index
        self.i7_index = i7_index
        self.i5_umi = re.findall(udelim, self.i5.sequence)
        self.i5_umi_before = 0
        self.i5_umi_after = 0
        self.i7_umi = re.findall(udelim, self.i7.sequence)
        self.i7_umi_before = 0
        self.i7_umi_after = 0
        self.name = f"{adaptor_type}"
        self.delim = delim

        if len(self.i5_umi) > 1 or len(self.i7_umi) > 1:
            raise UserWarning(
                f"Adaptor {adaptor_type} has more than one UMI in either i5 or i7. This is not supported."
            )
        # Check if UMI is before or after i5 index
        if len(self.i5_umi) > 0 and ">" + self.i5_umi[0] in self.i5.sequence:
            self.i5_umi_after = int(re.search(ulen, self.i5_umi[0]).group(1))
        elif len(self.i5_umi) > 0 and self.i5_umi[0] + "<" in self.i5.sequence:
            self.i5_umi_before = int(re.search(ulen, self.i5_umi[0]).group(1))
        elif len(self.i5_umi) > 0:
            raise UserWarning(
                f"Adaptor {adaptor_type} has UMI but it does not flank an index. This is not supported."
            )
        # Check if UMI is before or after i7 index
        if len(self.i7_umi) > 0 and ">" + self.i7_umi[0] in self.i7.sequence:
            self.i7_umi_after = int(re.search(ulen, self.i7_umi[0]).group(1))
        elif len(self.i7_umi) > 0 and self.i7_umi[0] + "<" in self.i7.sequence:
            self.i7_umi_before = int(re.search(ulen, self.i7_umi[0]).group(1))
        elif len(self.i7_umi) > 0:
            raise UserWarning(
                f"Adaptor {adaptor_type} has UMI but it does not flank an index. This is not supported."
            )

    def get_i5_mask(self, insert_Ns=True):
        ilen = len(self.i5_index) if self.i5_index is not None and insert_Ns else 0
        ulen = max(self.i5_umi_after, self.i5_umi_before) if insert_Ns else 0
        # Test if the index is specified in the adaptor sequence when it shouldn't be
        if has_match(idelim, self.i5.sequence) and self.i5_index is None and insert_Ns:
            raise UserWarning("Adaptor has i5 but no sequence was specified")
        if self.i5_index is not None or not insert_Ns:
            new_i5 = re.sub(idelim, "N" * ilen, self.i5.sequence)
            new_i5 = re.sub(udelim, "N" * ulen, new_i5)
            return new_i5
        else:
            return self.i5.sequence

    def get_i7_mask(self, insert_Ns=True):
        ilen = len(self.i7_index) if self.i7_index is not None and insert_Ns else 0
        ulen = max(self.i7_umi_after, self.i7_umi_before) if insert_Ns else 0
        # Test if the index is specified in the adaptor sequence when it shouldn't be
        if has_match(idelim, self.i7.sequence) and self.i7_index is None and insert_Ns:
            raise UserWarning("Adaptor has i7 but no sequence was specified")
        if self.i7_index is not None or not insert_Ns:
            new_i7 = re.sub(idelim, "N" * ilen, self.i7.sequence)
            new_i7 = re.sub(udelim, "N" * ulen, new_i7)
            return new_i7
        else:
            return self.i7.sequence

    def get_fastastring(self, insert_Ns=True):
        fasta_i5 = f">{self.name}_i5\n{self.get_i5_mask(insert_Ns)}\n"
        fasta_i7 = f">{self.name}_i7\n{self.get_i7_mask(insert_Ns)}\n"
        return fasta_i5 + fasta_i7


class AdaptorPart:
    # This class is used either the i5 or i7 adaptor
    def __init__(self, sequence, name, delim, index):
        self.sequence = sequence
        self.name = name
        self.delim = delim
        self.index = index
        self.umi_after = 0
        self.umi_before = 0
        self.len_after_index = 0
        self.len_before_index = 0

        # Dynamically assign attributes
        self.umi = re.findall(udelim, self.sequence)

        # TODO Duplicated from Adaptor class, will be merged later
        # Check if UMI is before or after index
        if len(self.umi) > 0 and ">" + self.umi[0] in self.sequence:
            # The index region is INDEX+UMI
            self.umi_after = int(re.search(ulen, self.umi[0]).group(1))
            self.len_before_index = len(idelim.split(self.sequence)[0])
            self.len_after_index = len(udelim.split(self.sequence)[-1])
        elif len(self.umi) > 0 and self.umi[0] + "<" in self.sequence:
            # The index region is UMI+INDEX
            self.umi_before = int(re.search(ulen, self.umi[0]).group(1))
            self.len_before_index = len(udelim.split(self.sequence)[0])
            self.len_after_index = len(idelim.split(self.sequence)[-1])
        elif len(self.umi) > 0:
            # TODO give details which adaptor has the problem
            raise UserWarning(
                "Found adaptor with UMI but it does not flank an index. This is not supported."
            )
        # Non UMI cases
        elif has_match(idelim, self.sequence):
            self.len_before_index = len(idelim.split(self.sequence)[0])
            self.len_after_index = len(idelim.split(self.sequence)[-1])

    def has_index(self):
        return self.sequence.find(self.delim) > -1

    def len_before_index_region(self):
        return self.len_before_index

    def len_after_index_region(self):
        return self.len_after_index


# General function to check if a string contains a pattern
def has_match(delim, seq):
    match = re.search(delim, seq)
    if match is None:
        return False
    return True


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
            Adaptor(adaptors_raw, "N", adaptor, i7_index=None, i5_index=None)
        )

    return adaptors
