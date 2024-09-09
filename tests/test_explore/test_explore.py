import gzip
import os
import tempfile

import pytest

from anglerfish.demux.adaptor import Adaptor
from anglerfish.explore.explore import (
    _run_explore,
    run_multiple_alignments,
    setup_explore_directory,
)


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

    # Read with perfect match to illumina ud i7 and i5 adaptors
    matching_read_seq = "".join(
        [
            "TTATGCTAACTCATTCAGCGTATTGCTAATGATACGGCGACCACCGAGATCTACAC",  # i5 before index
            "GCGTAAGA",  # i5 index
            "TCGTCGGCAGCGTC",  # i5 after index
            "CGGAAAATGTGGAAGTGACTTTGGAACTGGGTGCAGGAAAGACTGAAACAGTTTAGAGGGCGAAGAGAACGCGAAAAATATGGTAAGTTTGGAAAAGATCGGACGTCGTGAGAAAAGAGTGTAAAAGATCAGATGTGTAGATCTCA",  # Real read
            "CCGAGCCCACGAGAC",  # rev. comp i7 before index
            "TCCTGAGC",  # rev. comp i7 index
            "ATCTCGTATGCCGTCTTCTGCTTGAGCAATACGTT",  # rev. comp i7 after index
        ]
    )

    matching_read = (
        "\n".join(
            [
                "@matching_read",
                matching_read_seq,
                "+",
                "%" * len(matching_read_seq),
            ]
        )
        + "\n"
    )

    # Read with perfect match to illumina_ud, i5 only, slight mismatch i7
    i5_only_read_seq = "".join(
        [
            "TTATGCTAACTCATTCAGCGTATTGCTAATGATACGGCGACCACCGAGATCTACAC",  # i5 before index
            "GCGTAAGA",  # i5 index
            "TCGTCGGCAGCGTC",  # i5 after index
            "CGGAAAATGTGGAAGTGACTTTGGAACTGGGTGCAGGAAAGACTGAAACAGTTTAGAGGGCGAAGAGAACGCGAAAAATATGGTAAGTTTGGAAAAGATCGGACGTCGTGAGAAAAGAGTGTAAAAGATCAGATGTGTAGATCTCA",  # Real read
            "CAGACCCCACGAGAC",  # rev. comp i7 before index with 2 mismatches
            "TCCTGAGC",  # rev. comp i7 index
            "ATCTCGTATGCCGTCTTCTGCTTGAGCAATACGTT",  # rev. comp i7 after index
        ]
    )
    i5_only_read = (
        "\n".join(
            [
                "@half_matching_read",
                i5_only_read_seq,
                "+",
                "%" * len(i5_only_read_seq),
            ]
        )
        + "\n"
    )
    demux_fixture = request.getfixturevalue("demux_fixture")
    tmp_path = demux_fixture["tmp_path"]
    demux_fixture["explore_reads"] = tmp_path / "testdata_explore.fastq.gz"

    with gzip.open(tmp_path / "testdata_explore.fastq.gz", "ab") as f:
        f.write(matching_read.encode())
        f.write(i5_only_read.encode())

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


def test_setup_explore_directory():
    tmp_dir = tempfile.TemporaryDirectory()
    dir_to_be_created = os.path.join(tmp_dir.name, "test_dir")
    setup_explore_directory(dir_to_be_created, False)
    assert os.path.exists(dir_to_be_created)
    tmp_dir.cleanup()


def test_setup_explore_directory_use_existing():
    tmp_dir = tempfile.TemporaryDirectory()
    dir_to_be_created = os.path.join(tmp_dir.name, "test_dir")
    os.makedirs(dir_to_be_created)
    # This method should run exit(1) if the directory already exists
    with pytest.raises(SystemExit):
        setup_explore_directory(dir_to_be_created, False)

    # save modification time of the directory
    dir_time = os.path.getmtime(dir_to_be_created)
    setup_explore_directory(dir_to_be_created, True)
    # Check that the directory was not recreated
    assert dir_time == os.path.getmtime(dir_to_be_created)

    tmp_dir.cleanup()


def test_run_multiple_alignments(explore_fixture):
    # Test a basic run
    tmp_dir = tempfile.TemporaryDirectory()
    tmp_path = tmp_dir.name

    return_list = run_multiple_alignments(
        fastq=explore_fixture["testdata_single"],
        outdir=tmp_path,
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
        outdir=tmp_path,
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
        outdir=tmp_path,
        threads=1,
        use_existing=False,
        adaptors=explore_fixture["adaptors"],
        minimap_b=1,
    )

    # Check that the alignment file was rerun
    assert alignment_time != os.path.getmtime(alignment_file)

    tmp_dir.cleanup()


def test_run_explore_functional_test(explore_fixture):
    """Test overall function of the explore command."""
    tmp_dir = tempfile.TemporaryDirectory()
    tmp_path = os.path.join(tmp_dir.name, "outdir")

    # Running with 2 reads, one perfect match to illumina_ud, one partial match to illumina_ud
    results, adaptors_included, entries, umi_threshold, kmer_length, outdir = (
        _run_explore(
            fastq=explore_fixture["explore_reads"],
            outdir=tmp_path,
            threads=1,
            use_existing=False,
            good_hit_threshold=0.9,
            insert_thres_low=4,
            insert_thres_high=30,
            minimap_b=4,
            min_hits_per_adaptor=1,  # Not default
            umi_threshold=11,
            kmer_length=2,
        )
    )

    assert results

    # Make sure one read matches both and one read matches only i5
    assert list(results["included_adaptors"].keys()) == ["illumina_ud"]
    assert results["included_adaptors"]["illumina_ud"]["i5"]["nr_good_hits"] == 2
    assert results["included_adaptors"]["illumina_ud"]["i7"]["nr_good_hits"] == 1

    # Correct index length found
    assert results["included_adaptors"]["illumina_ud"]["i5"]["index_length"] == 8
    assert results["included_adaptors"]["illumina_ud"]["i7"]["index_length"] == 8

    assert len(adaptors_included) == 1
    assert adaptors_included[0].name == "illumina_ud"
