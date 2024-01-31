import importlib
import os
import re

import yaml

idelim = re.compile(r"\<N\>")
udelim = re.compile(r"(\<U\d+\>)")
ulen = re.compile(r"\<U(\d+)\>")


class Adaptor:
    def __init__(self, adaptors, delim, adaptor, i7_index=None, i5_index=None):
        self.i5 = AdaptorPart(adaptors[adaptor]["i7"], adaptor, delim, i7_index)
        self.i7 = AdaptorPart(adaptors[adaptor]["i5"], adaptor, delim, i5_index)
        self.i5_index = i5_index
        self.i7_index = i7_index
        self.i5_umi = re.findall(udelim, self.i5.sequence)
        self.i5_umi_before = 0
        self.i5_umi_after = 0
        self.i7_umi = re.findall(udelim, self.i7.sequence)
        self.i7_umi_before = 0
        self.i7_umi_after = 0
        self.name = f"{adaptor}"
        self.delim = delim

        if len(self.i5_umi) > 1 or len(self.i7_umi) > 1:
            raise UserWarning(
                f"Adaptor {adaptor} has more than one UMI in either i5 or i7. This is not supported."
            )
        # Check if UMI is before or after i5 index
        if len(self.i5_umi) > 0 and ">" + self.i5_umi[0] in self.i5.sequence:
            self.i5_umi_after = int(re.search(ulen, self.i5_umi[0]).group(1))
        elif len(self.i5_umi) > 0 and self.i5_umi[0] + "<" in self.i5.sequence:
            self.i5_umi_before = int(re.search(ulen, self.i5_umi[0]).group(1))
        elif len(self.i5_umi) > 0:
            raise UserWarning(
                f"Adaptor {adaptor} has UMI but it does not flank an index. This is not supported."
            )
        # Check if UMI is before or after i7 index
        if len(self.i7_umi) > 0 and ">" + self.i7_umi[0] in self.i7.sequence:
            self.i7_umi_after = int(re.search(ulen, self.i7_umi[0]).group(1))
        elif len(self.i7_umi) > 0 and self.i7_umi[0] + "<" in self.i7.sequence:
            self.i7_umi_before = int(re.search(ulen, self.i7_umi[0]).group(1))
        elif len(self.i7_umi) > 0:
            raise UserWarning(
                f"Adaptor {adaptor} has UMI but it does not flank an index. This is not supported."
            )
        if re.search(idelim, self.i5.sequence) is not None and i5_index is None:
            raise UserWarning("Adaptor has i5 but no sequence was specified")
        if re.search(idelim, self.i7.sequence) is not None and i7_index is None:
            raise UserWarning("Adaptor has i7 but no sequence was specified")

    def get_i5_mask(self):
        if self.i5_index is not None:
            new_i5 = re.sub(idelim, "N" * len(self.i5_index), self.i5.sequence)
            new_i5 = re.sub(
                udelim, "N" * max(self.i5_umi_after, self.i5_umi_before), new_i5
            )
            return new_i5
        else:
            return self.i5

    def get_i7_mask(self):
        if self.i7_index is not None:
            new_i7 = re.sub(idelim, "N" * len(self.i7_index), self.i7.sequence)
            new_i7 = re.sub(
                udelim, "N" * max(self.i7_umi_after, self.i7_umi_before), new_i7
            )
            return new_i7
        else:
            return self.i7


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
        # This is now broken, I think
        adaptors.append(
            Adaptor(adaptors_raw, "<N>", adaptor, i7_index=None, i5_index=None)
        )

    return adaptors
