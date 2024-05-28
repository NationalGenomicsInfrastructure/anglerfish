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
        self.index_token: str = INDEX_TOKEN
        self.i5 = AdaptorPart(
            sequence_token=adaptors[name]["i5"], name=name, end="i5", index_seq=i5_index
        )
        self.i7 = AdaptorPart(
            sequence_token=adaptors[name]["i7"], name=name, end="i7", index_seq=i7_index
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

    def __init__(self, sequence_token: str, name: str, end: str, index_seq: str | None):
        # Assign attributes from args
        self.sequence_token: str = sequence_token
        self.name: str = name
        self.end: str = end
        self.index_seq: str = index_seq

        # If index
        if has_match(INDEX_TOKEN, self.sequence_token):
            self.has_index = True
            self.len_index = len(index_seq) if index_seq else None
        else:
            self.has_index = False
            self.len_index = 0
            self.len_before_index = 0
            self.len_after_index = 0

        # UMI
        assert (
            len(re.findall(UMI_TOKEN, self.sequence_token)) <= 1
        ), f"Multiple UMIs found in {self.name}, not supported."
        if len(re.findall(UMI_TOKEN, self.sequence_token)) > 0:
            self.has_umi = True
            self.len_umi = int(
                re.search(UMI_LENGTH_TOKEN, self.sequence_token).group(1)
            )
        else:
            self.has_umi = False
            self.len_umi = 0
            self.len_umi_before_index = 0
            self.len_umi_after_index = 0

        if self.has_index and self.has_umi:
            # Check if UMI is before or after index
            if re.search(INDEX_TOKEN + UMI_TOKEN, self.sequence_token):
                # The index region is INDEX+UMI
                m = re.search(INDEX_TOKEN + UMI_TOKEN, self.sequence_token)
                self.len_umi_before_index = 0
                self.len_umi_after_index = self.len_umi
                self.len_before_index = len(self.sequence_token[: m.start()])
                self.len_after_index = (
                    len(self.sequence_token[m.end() :]) + self.len_umi
                )

            elif re.search(UMI_TOKEN + INDEX_TOKEN, self.sequence_token):
                # The index region is UMI+INDEX
                m = re.search(UMI_TOKEN + INDEX_TOKEN, self.sequence_token)
                self.len_umi_before_index = self.len_umi
                self.len_umi_after_index = 0
                self.len_before_index = (
                    len(self.sequence_token[: m.start()]) + self.len_umi
                )
                self.len_after_index = len(self.sequence_token[m.end() :])

            else:
                raise UserWarning(
                    f"Found adaptor {self.name} with UMI but it does not flank an index. This is not supported."
                )
        elif self.has_index and not self.has_umi:
            self.len_umi_before_index = 0
            self.len_umi_after_index = 0
            self.len_before_index = len(
                self.sequence_token[
                    : re.search(INDEX_TOKEN, self.sequence_token).start()
                ]
            )
            self.len_after_index = len(
                self.sequence_token[re.search(INDEX_TOKEN, self.sequence_token).end() :]
            )

    def get_mask(self, insert_Ns: bool = True) -> str:
        """Get the mask of the adaptor part.

        insert_Ns = True  -> Returns the sequence with index and UMI tokens replaced with Ns of the correct length
        insert_Ns = False -> Returns the sequence without index and UMI tokens
        """
        index_length = (
            len(self.index_seq) if self.index_seq is not None and insert_Ns else 0
        )
        umi_length = (
            max(self.len_umi_after_index, self.len_umi_before_index) if insert_Ns else 0
        )

        # Test if the index is specified in the adaptor sequence when it shouldn't be
        if (
            has_match(INDEX_TOKEN, self.sequence_token)
            and self.index_seq is None
            and insert_Ns
        ):
            raise UserWarning("Adaptor has i5 but no sequence was specified")

        if self.index_seq is not None or not insert_Ns:
            mask = re.sub(INDEX_TOKEN, "N" * index_length, self.sequence_token)
            mask = re.sub(UMI_TOKEN, "N" * umi_length, mask)
            return mask
        else:
            return self.sequence_token


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
