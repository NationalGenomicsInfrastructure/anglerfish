import re

from anglerfish.demux import adaptor as to_test


def test_has_match():
    assert (
        to_test.has_match(re.compile("pattern"), "thereisapatterninthisstring") is True
    )
    assert to_test.has_match(re.compile("pattern"), "nothinghere") is False


def test_validate_adaptors():
    valid_adaptors_dict = {
        "truseq_umi": {
            "i5": "TAGC<N>CCTT",
            "i7": "CAGT<N><U9>GGAA",
        },
    }
    assert to_test.validate_adaptors(valid_adaptors_dict) is None

    invalid_adaptors_dict = {
        "falseq_umi": {
            "i5": "TAGC<N>CCTT",
            "i7": "CAGT<N><U9>XYZ",
        },
    }
    try:
        to_test.validate_adaptors(invalid_adaptors_dict)
    except UserWarning as e:
        assert e
    else:
        raise AssertionError("UserWarning not raised")


def test_load_adaptors():
    adaptors = to_test.load_adaptors()
    assert adaptors
    assert isinstance(adaptors, list)
    for adaptor in adaptors:
        assert isinstance(adaptor, to_test.Adaptor)

    adaptors_raw = to_test.load_adaptors(raw=True)
    assert isinstance(adaptors_raw, dict)
    for adaptor in adaptors_raw:
        assert isinstance(adaptors_raw[adaptor], dict)


class TestAdaptorPart:
    """Explicit combinatorial testing, ugly but effective and readable.

    Here-in are contained test cases for a variety of instantiated AdaptorPart objects.
    All attributes and methods are tested for correctness.
    """

    def test_should_fail(self):
        """Specifying an index on an adaptor without an index should raise a UserWarning."""
        try:
            to_test.AdaptorPart(
                sequence_token="ATCG", name="should_fail", index_seq="AAA"
            )
        except UserWarning as e:
            assert e
        else:
            raise AssertionError("UserWarning not raised")

    def test_simple(self):
        adaptor_part = to_test.AdaptorPart(
            sequence_token="ATCG", name="simple", index_seq=None
        )
        assert adaptor_part.sequence_token == "ATCG"
        assert adaptor_part.name == "simple"
        assert adaptor_part.index_seq is None
        assert adaptor_part.has_index is False
        assert adaptor_part.len_index == 0
        assert adaptor_part.has_umi is False
        assert adaptor_part.len_umi == 0
        assert adaptor_part.len_before_index is None
        assert adaptor_part.len_after_index is None
        assert adaptor_part.len_umi_before_index is None
        assert adaptor_part.len_umi_after_index is None
        assert adaptor_part.len_total == 4
        assert adaptor_part.len_constant == 4
        assert adaptor_part.get_mask(insert_Ns=True) == "ATCG"
        assert adaptor_part.get_mask(insert_Ns=False) == "ATCG"

    def test_unknown_index(self):
        adaptor_part = to_test.AdaptorPart(
            sequence_token="ATCG<N>ATC", name="unknown_index", index_seq=None
        )
        assert adaptor_part.sequence_token == "ATCG<N>ATC"
        assert adaptor_part.name == "unknown_index"
        assert adaptor_part.index_seq is None
        assert adaptor_part.has_index is True
        assert adaptor_part.len_index is None
        assert adaptor_part.has_umi is False
        assert adaptor_part.len_umi == 0
        assert adaptor_part.len_before_index == 4
        assert adaptor_part.len_after_index == 3
        assert adaptor_part.len_umi_before_index == 0
        assert adaptor_part.len_umi_after_index == 0
        assert adaptor_part.len_total is None
        assert adaptor_part.len_constant == 7
        try:
            adaptor_part.get_mask(insert_Ns=True)
        except UserWarning as e:
            assert e
        else:
            raise AssertionError("UserWarning not raised")
        assert adaptor_part.get_mask(insert_Ns=False) == "ATCGATC"

    def test_known_index(self):
        adaptor_part = to_test.AdaptorPart(
            sequence_token="ATCG<N>ATC", name="known_index", index_seq="GGG"
        )
        assert adaptor_part.sequence_token == "ATCG<N>ATC"
        assert adaptor_part.name == "known_index"
        assert adaptor_part.index_seq == "GGG"
        assert adaptor_part.has_index is True
        assert adaptor_part.len_index == 3
        assert adaptor_part.has_umi is False
        assert adaptor_part.len_umi == 0
        assert adaptor_part.len_before_index == 4
        assert adaptor_part.len_after_index == 3
        assert adaptor_part.len_umi_before_index == 0
        assert adaptor_part.len_umi_after_index == 0
        assert adaptor_part.len_total == 10
        assert adaptor_part.len_constant == 7
        assert adaptor_part.get_mask(insert_Ns=True) == "ATCGNNNATC"
        assert adaptor_part.get_mask(insert_Ns=False) == "ATCGATC"

    def test_unknown_index_umi(self):
        adaptor_part = to_test.AdaptorPart(
            sequence_token="ATCG<N><U5>ATC", name="unknown_index_umi", index_seq=None
        )
        assert adaptor_part.sequence_token == "ATCG<N><U5>ATC"
        assert adaptor_part.name == "unknown_index_umi"
        assert adaptor_part.index_seq is None
        assert adaptor_part.has_index is True
        assert adaptor_part.len_index is None
        assert adaptor_part.has_umi is True
        assert adaptor_part.len_umi == 5
        assert adaptor_part.len_before_index == 4
        assert adaptor_part.len_after_index == 8
        assert adaptor_part.len_umi_before_index == 0
        assert adaptor_part.len_umi_after_index == 5
        assert adaptor_part.len_total is None
        assert adaptor_part.len_constant == 7
        try:
            adaptor_part.get_mask(insert_Ns=True)
        except UserWarning as e:
            assert e
        else:
            raise AssertionError("UserWarning not raised")
        assert adaptor_part.get_mask(insert_Ns=False) == "ATCGATC"

    def test_known_index_umi(self):
        adaptor_part = to_test.AdaptorPart(
            sequence_token="ATCG<N><U5>ATC", name="known_index_umi", index_seq="GGG"
        )
        assert adaptor_part.sequence_token == "ATCG<N><U5>ATC"
        assert adaptor_part.name == "known_index_umi"
        assert adaptor_part.index_seq == "GGG"
        assert adaptor_part.has_index is True
        assert adaptor_part.len_index == 3
        assert adaptor_part.has_umi is True
        assert adaptor_part.len_umi == 5
        assert adaptor_part.len_before_index == 4
        assert adaptor_part.len_after_index == 8
        assert adaptor_part.len_umi_before_index == 0
        assert adaptor_part.len_umi_after_index == 5
        assert adaptor_part.len_total == 15
        assert adaptor_part.len_constant == 7
        assert adaptor_part.get_mask(insert_Ns=True) == "ATCGNNNNNNNNATC"
        assert adaptor_part.get_mask(insert_Ns=False) == "ATCGATC"

    def test_umi_known_index(self):
        adaptor_part = to_test.AdaptorPart(
            sequence_token="ATCG<U5><N>ATC", name="umi_known_index", index_seq="GGG"
        )
        assert adaptor_part.sequence_token == "ATCG<U5><N>ATC"
        assert adaptor_part.name == "umi_known_index"
        assert adaptor_part.index_seq == "GGG"
        assert adaptor_part.has_index is True
        assert adaptor_part.len_index == 3
        assert adaptor_part.has_umi is True
        assert adaptor_part.len_umi == 5
        assert adaptor_part.len_before_index == 9
        assert adaptor_part.len_after_index == 3
        assert adaptor_part.len_umi_before_index == 5
        assert adaptor_part.len_umi_after_index == 0
        assert adaptor_part.len_total == 15
        assert adaptor_part.len_constant == 7
        assert adaptor_part.get_mask(insert_Ns=True) == "ATCGNNNNNNNNATC"
        assert adaptor_part.get_mask(insert_Ns=False) == "ATCGATC"


class TestAdaptor:
    def test_adaptor(self):
        adaptors = {
            "simple_and_index_umi": {
                "i5": "AAA",
                "i7": "AAA<N><U4>CCC",
            }
        }

        adaptor = to_test.Adaptor(
            "simple_and_index_umi", adaptors, i5_index=None, i7_index="TTT"
        )
        assert adaptor.name == "simple_and_index_umi"
        assert isinstance(adaptor.i5, to_test.AdaptorPart)
        assert isinstance(adaptor.i7, to_test.AdaptorPart)
        assert (
            adaptor.get_fastastring(insert_Ns=True)
            == "\n".join(
                [
                    ">simple_and_index_umi_i5",
                    "AAA",
                    ">simple_and_index_umi_i7",
                    "AAANNNNNNNCCC",
                ]
            )
            + "\n"
        )
        assert (
            adaptor.get_fastastring(insert_Ns=False)
            == "\n".join(
                [
                    ">simple_and_index_umi_i5",
                    "AAA",
                    ">simple_and_index_umi_i7",
                    "AAACCC",
                ]
            )
            + "\n"
        )
