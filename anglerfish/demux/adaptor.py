import importlib.resources as resources
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
        self.i5 = AdaptorPart(
            sequence_token=adaptors[name]["i5"], name=name, index_seq=i5_index
        )
        self.i7 = AdaptorPart(
            sequence_token=adaptors[name]["i7"], name=name, index_seq=i7_index
        )

    def get_fastastring(self, insert_Ns: bool = True) -> str:
        """Create fasta text entries for the i5 and i7 adaptors.

        Example:

            An Adaptor object with the name "truseq_dual" and the following properties:
                - i5 sequence token 'AAAA<N>TTTT', index length 4
                - i7 sequence token 'GGGG<N><U8>CCCC', index length 4, UMI length 8

            get_fastastring(insert_Ns=True) -> '''
                >truseq_dual_i5
                AAAANNNNTTTT
                >truseq_dual_i7
                GGGGNNNNNNNNNNNNCCCC
            '''

            get_fastastring(insert_Ns=False) -> '''
                >truseq_dual_i5
                AAAATTTT
                >truseq_dual_i7
                GGGGCCCC
            '''

        """
        fasta_i5 = f">{self.name}_i5\n{self.i5.get_mask(insert_Ns)}\n"
        fasta_i7 = f">{self.name}_i7\n{self.i7.get_mask(insert_Ns)}\n"
        return fasta_i5 + fasta_i7


class AdaptorPart:
    """This class is used for the i5 or i7 adaptor."""

    def __init__(self, sequence_token: str, name: str, index_seq: str | None):
        ## Type declaration of attributes to be assigned upon instantiation
        # Attributes from arguments
        self.name: str
        self.sequence_token: str
        self.index_seq: str | None

        # Index attributes
        self.has_index: bool
        self.len_index: int | None
        self.len_before_index: int | None
        self.len_after_index: int | None

        # UMI attributes
        self.has_umi: bool
        self.len_umi: int | None
        self.len_umi_before_index: int | None
        self.len_umi_after_index: int | None

        # Length attributes
        self.len_total: int | None
        self.len_constant: int

        # Instantiation outsorced to private method
        self._setup(sequence_token, name, index_seq)

    def _setup(self, sequence_token: str, name: str, index_seq: str | None):
        # Assign attributes from args
        self.sequence_token = sequence_token
        self.name = name
        self.index_seq = index_seq

        # Index bool and len
        if has_match(INDEX_TOKEN, self.sequence_token):
            split_by_index = re.split(INDEX_TOKEN, self.sequence_token)

            self.has_index = True
            self.len_index = len(index_seq) if index_seq else None

        else:
            if self.index_seq is not None:
                raise UserWarning(
                    "Index sequence specified, but no index token found in adaptor sequence."
                )
            self.has_index = False
            self.len_index = 0

        # UMI bool and len
        umi_tokens = re.findall(UMI_TOKEN, self.sequence_token)
        if len(umi_tokens) > 1:
            raise UserWarning(
                f"Found adaptor {self.name} with multiple UMIs. This is not supported."
            )
        elif len(umi_tokens) == 1:
            self.has_umi = True
            umi_token_search = re.search(UMI_LENGTH_TOKEN, self.sequence_token)
            assert isinstance(umi_token_search, re.Match)
            self.len_umi = int(umi_token_search.group(1))
        else:
            self.has_umi = False
            self.len_umi = 0

        # Lengths
        if self.has_index and self.has_umi:
            # Index and UMI
            index_umi_match = re.search(INDEX_TOKEN + UMI_TOKEN, self.sequence_token)
            umi_index_match = re.search(UMI_TOKEN + INDEX_TOKEN, self.sequence_token)

            if index_umi_match:
                self.len_umi_before_index = 0
                self.len_umi_after_index = self.len_umi
                self.len_before_index = len(
                    self.sequence_token[: index_umi_match.start()]
                )
                self.len_after_index = (
                    len(self.sequence_token[index_umi_match.end() :]) + self.len_umi
                )
            elif umi_index_match:
                self.len_umi_before_index = self.len_umi
                self.len_umi_after_index = 0
                self.len_before_index = (
                    len(self.sequence_token[: umi_index_match.start()]) + self.len_umi
                )
                self.len_after_index = len(self.sequence_token[umi_index_match.end() :])
            else:
                raise UserWarning(
                    f"Found adaptor {self.name} with UMI but it does not flank an index. This is not supported."
                )

        elif self.has_index and not self.has_umi:
            # Index, no UMI
            self.len_umi_before_index = 0
            self.len_umi_after_index = 0
            self.len_before_index = len(split_by_index[0])
            self.len_after_index = len(split_by_index[-1])

        elif not self.has_index and self.has_umi:
            # No index, UMI
            raise UserWarning(
                f"Adaptor {self.name} has UMI but no index. This is not supported."
            )

        else:
            # No index, no UMI
            self.len_umi_before_index = None
            self.len_umi_after_index = None
            self.len_before_index = None
            self.len_after_index = None

        if (
            self.has_index is True and self.index_seq is not None
        ) or self.has_index is False:
            self.len_total = len(self.get_mask(insert_Ns=True))
        else:
            self.len_total = None
        self.len_constant = len(self.get_mask(insert_Ns=False))

    def get_mask(self, insert_Ns: bool = True) -> str:
        """Get the mask of the adaptor part.

        insert_Ns = True  -> Returns the sequence with index and UMI tokens replaced with Ns of the correct length
        insert_Ns = False -> Returns the sequence without index and UMI tokens
        """

        index_mask_length = (
            len(self.index_seq)
            if self.index_seq is not None and insert_Ns and self.has_index
            else 0
        )

        if insert_Ns and self.has_umi:
            assert self.len_umi_before_index is not None
            assert self.len_umi_after_index is not None
            umi_mask_length = max(self.len_umi_after_index, self.len_umi_before_index)
        else:
            umi_mask_length = 0

        # Test if the index is specified in the adaptor sequence when it shouldn't be
        if (
            has_match(INDEX_TOKEN, self.sequence_token)
            and self.index_seq is None
            and insert_Ns
        ):
            raise UserWarning(
                f"Can't create mask for adaptor '{self.name}' with unspecified index."
            )

        if self.index_seq is not None or not insert_Ns:
            mask = re.sub(INDEX_TOKEN, "N" * index_mask_length, self.sequence_token)
            mask = re.sub(UMI_TOKEN, "N" * umi_mask_length, mask)
            return mask
        else:
            return self.sequence_token


def has_match(pattern: re.Pattern | str, query: str) -> bool:
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

    return True


def load_adaptors(raw: bool = False) -> list[Adaptor] | dict:
    """Fetch all adaptors.

    Return them as a list of adaptor objects or optionally as a raw yaml dict.
    """

    # Load adaptors from config file
    adaptors_config_path = resources.files("anglerfish.config").joinpath(
        "adaptors.yaml"
    )
    assert isinstance(adaptors_config_path, os.PathLike)

    with open(adaptors_config_path) as f:
        adaptors_dict = yaml.safe_load(f)

    # Validate input
    assert validate_adaptors(adaptors_dict) is True

    # Optionally, return raw dict
    if raw:
        return adaptors_dict

    # By default, return list of Adaptor objects
    else:
        adaptors = []
        for adaptor_name in adaptors_dict:
            adaptors.append(Adaptor(name=adaptor_name, adaptors=adaptors_dict))
        return adaptors
