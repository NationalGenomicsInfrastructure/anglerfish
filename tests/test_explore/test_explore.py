import os

import pytest

from anglerfish.demux.adaptor import Adaptor
from anglerfish.explore.explore import run_multiple_alignments


@pytest.fixture(scope="module")
def explore_fixture(request):
    f"""Fixture for all tests within {__file__}.
    
    Creates a temporary directory with...
    - truseq.fasta: A fasta file with the truseq i7 and i5 adaptors.
    - illumina_ud.fasta: A fasta file with the illumina unique dual i7 and i5 adaptors.
    - testdata_single.fastq.gz: A single read subset from the testdata with a perfect index match (TAACTTGGTC) to the truseq i7 adaptor.
    - testdata_multiple.fastq.gz: Fabricated test data intended to test the categorization of alignments.

    Expected files to be created...
    """

    adaptors_dict = {
        "illumina_ud": {
            "i5": "AATGATACGGCGACCACCGAGATCTACAC<N>TCGTCGGCAGCGTC",
            "i7": "CAAGCAGAAGACGGCATACGAGAT<N>GTCTCGTGGGCTCGG",
        },
        "truseq": {
            "i5": "AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT",
            "i7": "GATCGGAAGAGCACACGTCTGAACTCCAGTCAC<N>ATCTCGTATGCCGTCTTCTGCTTG",
        },
    }
    demux_fixture = request.getfixturevalue("demux_fixture")
    fixture = demux_fixture
    fixture["adaptors"] = [
        Adaptor(
            name=adaptor_name,
            i7_sequence_token=adaptor_parts_dict["i7"],
            i5_sequence_token=adaptor_parts_dict["i5"],
        )
        for adaptor_name, adaptor_parts_dict in adaptors_dict.items()
    ]
    return fixture


def test_run_multiple_alignments(explore_fixture):
    # Test a basic run
    return_list = run_multiple_alignments(
        fastq=explore_fixture["testdata_single"],
        outdir=explore_fixture["tmp_path"],
        threads=1,
        use_existing=False,
        adaptors=explore_fixture["adaptors"],
        minimap_b=1,
    )

    assert return_list
    assert len(return_list) == 2  # One for illumina_ud and one for truseq
    for adapter, alignment_file in return_list:
        # check that tile file ending is .paf
        assert alignment_file.endswith(".paf")
        # check that alignment file exists
        assert os.path.exists(alignment_file)

    # When use_existing is true, the alignment should not be rerun
    # Save modification time of the alignment file
    alignment_file = return_list[0][1]
    alignment_time = os.path.getmtime(alignment_file)

    return_list = run_multiple_alignments(
        fastq=explore_fixture["testdata_single"],
        outdir=explore_fixture["tmp_path"],
        threads=1,
        use_existing=True,
        adaptors=explore_fixture["adaptors"],
        minimap_b=1,
    )
    # Check that the alignment file was not rerun
    assert alignment_time == os.path.getmtime(alignment_file)

    # When use_existing is false, the alignment should be rerun
    # Save modification time of the alignment file
    alignment_file = return_list[0][1]
    alignment_time = os.path.getmtime(alignment_file)

    return_list = run_multiple_alignments(
        fastq=explore_fixture["testdata_single"],
        outdir=explore_fixture["tmp_path"],
        threads=1,
        use_existing=False,
        adaptors=explore_fixture["adaptors"],
        minimap_b=1,
    )
    # Check that the alignment file was rerun
    assert alignment_time != os.path.getmtime(alignment_file)
