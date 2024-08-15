import gzip
import tempfile
from pathlib import Path

import pytest

from anglerfish.demux import demux as to_test


@pytest.fixture(scope="module")
def fixture():
    f"""Fixture for all tests within {__file__}.
    
    Creates a temporary directory with...
    - truseq.fasta: A fasta file with the truseq i7 and i5 adaptors.
    - testdata_single.fastq.gz: A single read subset from the testdata with a perfect index match (TAACTTGGTC) to the truseq i7 adaptor.
    - testdata_multiple.fastq.gz: Fabricated test data intended to test the categorization of alignments.

    Expected files to be created...
    - test_single.paf: The output of minimap2 run on testdata_single.fastq.gz and truseq.fasta.
    - test_multiple.paf: The output of minimap2 run on testdata_multiple.fastq.gz and truseq.fasta.
    """
    tmp = tempfile.TemporaryDirectory()
    tmp_path = Path(tmp.name)

    # Create testdata of single read
    with gzip.open(tmp_path / "testdata_single.fastq.gz", "wb") as f:
        contents = "\n".join(
            [
                "@0ad8bdb6-e009-43c5-95b1-d381e699f983 runid=7d98ab9ee7f513f92fd5e382dc3dd563b956e4fa sampleid=Pools read=7283 ch=21 start_time=2020-01-17T15:19:22Z",
                "CGTTGTAGTTCACCAAACCCAACAACCTAGATAGGCCGCACCTAAAATGATACGGCGACCGCGAGATCCCGATTCTGACACTCTTTCCCTACACGACGCTCTTCCGATCTCACAAAATTATTAGGCATGGTGGTGCATGCCTGTGGTCCCAGCTACTTGGGAGGCTGAGGCAGGGAATCGCTTGAACCCAGGAGGCAGAGGTTGCAGTGAGCCAGACCAGCGCCACTGTACTCCAGCCTGGCAACAAGCAGGAACTCCGTTTCAAAACAAGATCTAGGGTGGTTTTCCCAAGTTAGTCAGGAGATCGGAAGAGCACGTCGGAACTCCAGTCACTAACTTGGTCAACCTAGAAATCTGATGCCGTCTTCTGCTTGTTAGGTGCTGGCCTATCTGGGTGTTGGGTTTGGTTGTATAACGT",
                "+",
                "&&##$$'$()'$$&$&%%%%%%%(*340++1','',,)-*669:98731/0+&&(%&$$$%(',,4$$,&;>@:2$-')&&25<<>/4.,'*$,--**6<;<>)1...2-&%+///+++*%'*1B:<BFDDB=760998EEA?%33115/5-&&3%%/22:+07660*76..7221680('+058)777(@.A?5.-+&'*7)1*-$)))-/1)/'$*%(**.6*),$&&,+%13-2+'%$(%&+&(&'+%$$+.0023<<6946>93;:8618'/&%,%'()344:=:'*6*5$*39090-/4/4)537-/166+)%,'(%'02)5-5@A:/;75&'3?=9=--=;.,/030'9%&&((48977;<=@;3445%'.2<=;51/(.,),*7*%&(,(()&'/01*)'%%$%%(%%$%$",
                "",
            ]
        )
        f.write(contents.encode())

    # Create testdata of multiple different reads
    dummy_read = "ATCG" * 10
    truseq_i5 = "AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT"
    truseq_index = "TAACTTGGTC"
    truseq_i7 = (
        f"GATCGGAAGAGCACACGTCTGAACTCCAGTCAC{truseq_index}ATCTCGTATGCCGTCTTCTGCTTG"
    )

    test_reads = {
        "fragment_read": truseq_i5 + dummy_read + truseq_i7,
        "singleton_i7_read": dummy_read + truseq_i7,
        "singleton_i5_read": truseq_i5 + dummy_read,
        "concat_read": truseq_i5 + truseq_i5 + dummy_read + truseq_i7 + truseq_i7,
        "unknown_read": truseq_i5 + dummy_read + truseq_i5,
    }

    with gzip.open(tmp_path / "testdata_multiple.fastq.gz", "ab") as f:
        for read_name, read_sequence in test_reads.items():
            contents = "\n".join(
                [f"@{read_name}", read_sequence, "+", len(read_sequence) * "%", ""]
            )
            f.write(contents.encode())

    # Create adaptor fasta file
    with open(tmp_path / "truseq.fasta", "w") as f:
        contents = "\n".join(
            [
                ">truseq_i7",
                "GATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNNNNNATCTCGTATGCCGTCTTCTGCTTG",
                ">truseq_i5",
                "AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT",
                "",
            ]
        )
        f.write(contents)

    fixture = {
        "tmp_path": tmp_path,
        "index_file": tmp_path / "truseq.fasta",
        "index": truseq_index,
        "testdata_single": tmp_path / "testdata_single.fastq.gz",
        "testdata_multiple": tmp_path / "testdata_multiple.fastq.gz",
        "paf_single": tmp_path / "test_single.paf",
        "paf_multiple": tmp_path / "test_multiple.paf",
    }

    yield fixture

    tmp.cleanup()


def test_run_minimap2(fixture):
    """Check that the function runs successfully, not the output."""
    # Test alignment on single read
    to_test.run_minimap2(
        fastq_in=fixture["testdata_single"],
        index_file=fixture["index_file"],
        output_paf=fixture["paf_single"],
        threads=1,
        minimap_b=1,
    )

    # Create aligntment from multiple reads
    to_test.run_minimap2(
        fastq_in=fixture["testdata_multiple"],
        index_file=fixture["index_file"],
        output_paf=fixture["paf_multiple"],
        threads=1,
        minimap_b=1,
    )


def test_map_reads_to_alns(fixture):
    reads_alns = to_test.map_reads_to_alns(fixture["paf_single"])

    for read_name, alns in reads_alns.items():
        assert read_name == "0ad8bdb6-e009-43c5-95b1-d381e699f983"
        for aln in alns:
            if aln.adapter_name == "truseq_i7":
                assert (
                    aln.cs
                    == "cs:Z::11-ca:6*tg:13*nt*na*na*nc*nt*nt*ng*ng*nt*nc:1*ta:1+ctagaaa:2*gt*tg:17"
                )
            else:
                assert aln.cs == "cs:Z::15-a*cg:5+tcccgat:3+g:33"


def test_parse_cs(fixture):
    test_cs_str = (
        "cs:Z::11-ca:6*tg:13*nt*na*na*nc*nt*nt*ng*ng*nt*nc:1*ta:1+ctagaaa:2*gt*tg:17"
    )
    expected_cs_parsed = ("taacttggtc", 0)

    cs_parsed = to_test.parse_cs(
        cs_string=test_cs_str,
        index_seq=fixture["index"],
        umi_before=0,
        umi_after=0,
    )

    assert cs_parsed == expected_cs_parsed


def test_categorize_matches(fixture):
    i5_name = "truseq_i5"
    i7_name = "truseq_i7"
    reads_alns = to_test.map_reads_to_alns(fixture["paf_multiple"])

    layout = to_test.categorize_matches(
        i5_name=i5_name, i7_name=i7_name, reads_to_alns=reads_alns
    )
    fragments, singletons, concats, unknowns = layout

    # Assert reads were categorized as expected
    assert len(fragments["fragment_read"]) == 2
    assert len(singletons["singleton_i5_read"]) == 1
    assert len(singletons["singleton_i7_read"]) == 1
    assert len(concats["concat_read"]) == 4
    assert len(unknowns["unknown_read"]) == 2
