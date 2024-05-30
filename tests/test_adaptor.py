import re

import pytest

from anglerfish.demux import adaptor as adaptor_to_test


def test_has_match():
    assert (
        adaptor_to_test.has_match(re.compile("pattern"), "thereisapatterninthisstring") is True
    )
    assert adaptor_to_test.has_match(re.compile("pattern"), "nothinghere") is False


def test_validate_adaptors():
    valid_adaptors_dict = {
        "truseq_umi": {
            "i5": "TAGC<N>CCTT",
            "i7": "CAGT<N><U9>GGAA",
        },
    }
    assert adaptor_to_test.validate_adaptors(valid_adaptors_dict) is None

    invalid_adaptors_dict = {
        "falseq_umi": {
            "i5": "TAGC<N>CCTT",
            "i7": "CAGT<N><U9>XYZ",
        },
    }
    try:
        adaptor_to_test.validate_adaptors(invalid_adaptors_dict)
    except UserWarning as e:
        assert e
    else:
        raise AssertionError("UserWarning not raised")


def test_load_adaptors():
    adaptors = adaptor_to_test.load_adaptors()
    assert adaptors
    assert isinstance(adaptors, list)
    for adaptor in adaptors:
        assert isinstance(adaptor, adaptor_to_test.Adaptor)
