import importlib
import os
import re

import yaml

# These variables correspond to the tokens used in the adaptors.yaml file
# Not compiled with re.compile to enable pattern concatenation
INDEX_TOKEN = r"(<N>)"
UMI_TOKEN = r"(\<U\d+\>)"
UMI_LENGTH_TOKEN = r"\<U(\d+)\>"

# This pattern is used to validate the adaptor sequences in the config
VALID_SEQUENCE_TOKEN_PATTERN = re.compile(f"^({INDEX_TOKEN}|{UMI_TOKEN}|([ACTG]*))*$")


class Adaptor:
    def __init__(
        self,
        name: str,
        adaptors: dict,
        i7_index: str | None = None,
        i5_index: str | None = None,
    ):
        self.name: str = name
        self.i5_token = (adaptors[name]["i5"],)
        self.i7_token = (adaptors[name]["i7"],)
        self.index_token: str = INDEX_TOKEN

        # i5 attributes
        self.i5 = AdaptorPart(
            sequence_token=self.i5_token,
            name=name,
            index=i5_index,
        )
        self.i5_index: str | None = i5_index
        self.i5_umi: str | None = self.i5.umi_token
        self.i5_umi_before: int = self.i5.len_umi_before_index
        self.i5_umi_after: int = self.i5.len_umi_after_index

        # i7 attributes
        self.i7 = AdaptorPart(
            sequence_token=self.i7_token,
            name=name,
            index=i7_index,
        )
        self.i7_index: str | None = i7_index
        self.i7_umi: str | None = self.i7.umi_token
        self.i7_umi_before: int = self.i7.len_umi_before_index
        self.i7_umi_after: int = self.i7.len_umi_after_index

    def get_i5_mask(self, insert_Ns: bool = True) -> str:
        ilen = len(self.i5_index) if self.i5_index is not None and insert_Ns else 0
        ulen = max(self.i5_umi_after, self.i5_umi_before) if insert_Ns else 0
        # Test if the index is specified in the adaptor sequence when it shouldn't be
        if (
            has_match(INDEX_TOKEN, self.i5.sequence_token)
            and self.i5_index is None
            and insert_Ns
        ):
            raise UserWarning("Adaptor has i5 but no sequence was specified")
        if self.i5_index is not None or not insert_Ns:
            new_i5 = re.sub(INDEX_TOKEN, "N" * ilen, self.i5.sequence_token)
            new_i5 = re.sub(UMI_TOKEN, "N" * ulen, new_i5)
            return new_i5
        else:
            return self.i5.sequence_token

    def get_i7_mask(self, insert_Ns: bool = True) -> str:
        ilen = len(self.i7_index) if self.i7_index is not None and insert_Ns else 0
        ulen = max(self.i7_umi_after, self.i7_umi_before) if insert_Ns else 0
        # Test if the index is specified in the adaptor sequence when it shouldn't be
        if (
            has_match(INDEX_TOKEN, self.i7.sequence_token)
            and self.i7_index is None
            and insert_Ns
        ):
            raise UserWarning("Adaptor has i7 but no sequence was specified")
        if self.i7_index is not None or not insert_Ns:
            new_i7 = re.sub(INDEX_TOKEN, "N" * ilen, self.i7.sequence_token)
            new_i7 = re.sub(UMI_TOKEN, "N" * ulen, new_i7)
            return new_i7
        else:
            return self.i7.sequence_token

    def get_fastastring(self, insert_Ns: bool = True) -> str:
        fasta_i5 = f">{self.name}_i5\n{self.get_i5_mask(insert_Ns)}\n"
        fasta_i7 = f">{self.name}_i7\n{self.get_i7_mask(insert_Ns)}\n"
        return fasta_i5 + fasta_i7


class AdaptorPart:
    """This class is used for the i5 or i7 adaptor."""

    def __init__(self, sequence_token: str, name: str, index: str | None):
        # Assign attributes from args
        self.sequence_token: str = sequence_token
        self.name: str = name
        self.index: str = index

        # Parse index, if any
        self.has_index: bool = has_match(INDEX_TOKEN, self.sequence_token)

        # Parse UMI, if any
        umi_token_matches = re.findall(UMI_TOKEN, self.sequence_token)
        if umi_token_matches > 0:
            assert (
                umi_token_matches == 1
            ), f"Multiple UMIs found in {self.name}, not supported."
            self.umi_token = umi_token_matches[0]

            # Check if UMI is before or after index
            if INDEX_TOKEN + UMI_TOKEN in self.sequence_token:
                # The index region is INDEX+UMI
                self.len_umi_after_index = int(
                    re.search(UMI_LENGTH_TOKEN, self.umi_token).group(1)
                )
                self.len_before_index = len(INDEX_TOKEN.split(self.sequence_token)[0])
                self.len_after_index = len(UMI_TOKEN.split(self.sequence_token)[-1])

            elif UMI_TOKEN + INDEX_TOKEN in self.sequence_token:
                # The index region is UMI+INDEX
                self.len_umi_before_index = int(
                    re.search(UMI_LENGTH_TOKEN, self.umi_token[0]).group(1)
                )
                self.len_before_index = len(UMI_TOKEN.split(self.sequence_token)[0])
                self.len_after_index = len(INDEX_TOKEN.split(self.sequence_token)[-1])

            else:
                raise UserWarning(
                    f"Found adaptor {self.name} with UMI but it does not flank an index. This is not supported."
                )

        else:
            self.umi_token = None
            self.len_before_index = len(INDEX_TOKEN.split(self.sequence_token)[0])
            self.len_after_index = len(INDEX_TOKEN.split(self.sequence_token)[-1])


def has_match(pattern: re.Pattern, query: str) -> bool:
    """General function to check if a string contains a pattern."""
    match = re.search(pattern, query)
    if match is None:
        return False
    return True


def validate_adaptors(adaptors_dict: dict):
    """Validate that the adaptor config sequences only consist of expected patterns."""

    for adaptor_name in adaptors_dict:
        for i in ["i5", "i7"]:
            sequence_token = adaptors_dict[adaptor_name][i]
            match = re.match(VALID_SEQUENCE_TOKEN_PATTERN, sequence_token)
            if not match:
                raise UserWarning(
                    f"Adaptor {adaptor_name} has an invalid sequence for {i}: {sequence_token}. Does not conform to the pattern {VALID_SEQUENCE_TOKEN_PATTERN}."
                )


def load_adaptors(raw: bool = False) -> list[Adaptor] | dict:
    """Fetch all adaptors.

    Return them as a list of adaptor objects or optionally as a raw yaml dict.
    """

    # Load adaptors from config file
    adaptors_config_path = importlib.resources.files("anglerfish.config").joinpath(
        "adaptors.yaml"
    )
    assert isinstance(adaptors_config_path, os.PathLike)

    with open(adaptors_config_path) as f:
        adaptors_dict = yaml.safe_load(f)

    # Validate input
    validate_adaptors(adaptors_dict)

    # Optionally, return raw dict
    if raw:
        return adaptors_dict

    # By default, return list of Adaptor objects
    else:
        adaptors = []
        for adaptor_name in adaptors_dict:
            adaptors.append(Adaptor(name=adaptor_name, adaptors=adaptors_dict))
        return adaptors
