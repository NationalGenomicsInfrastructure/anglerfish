import gzip
import tempfile
from pathlib import Path

import pytest

from anglerfish.demux import demux as to_test


@pytest.fixture(scope="module")
def fixture():
    f"""Fixture for all tests within {__file__}.
    
    Creats a temporary directory with...
    - testdata.fastq.gz: A single read subset from the testdata with a perfect index match to the truseq i7 adaptor.
    - truseq.fasta: A fasta file with the truseq i7 and i5 adaptors.


    Expected files to be created...
    - test.paf: The output of minimap2 run on the testdata.fastq.gz and truseq.fasta.
    """
    tmp = tempfile.TemporaryDirectory()
    tmp_path = Path(tmp.name)

    # Create test reads
    with gzip.open(tmp_path / "testdata.fastq.gz", "wb") as f:
        f.write(b"""
@0ad8bdb6-e009-43c5-95b1-d381e699f983 runid=7d98ab9ee7f513f92fd5e382dc3dd563b956e4fa sampleid=Pools read=7283 ch=21 start_time=2020-01-17T15:19:22Z
CGTTGTAGTTCACCAAACCCAACAACCTAGATAGGCCGCACCTAAAATGATACGGCGACCGCGAGATCCCGATTCTGACACTCTTTCCCTACACGACGCTCTTCCGATCTCACAAAATTATTAGGCATGGTGGTGCATGCCTGTGGTCCCAGCTACTTGGGAGGCTGAGGCAGGGAATCGCTTGAACCCAGGAGGCAGAGGTTGCAGTGAGCCAGACCAGCGCCACTGTACTCCAGCCTGGCAACAAGCAGGAACTCCGTTTCAAAACAAGATCTAGGGTGGTTTTCCCAAGTTAGTCAGGAGATCGGAAGAGCACGTCGGAACTCCAGTCACTAACTTGGTCAACCTAGAAATCTGATGCCGTCTTCTGCTTGTTAGGTGCTGGCCTATCTGGGTGTTGGGTTTGGTTGTATAACGT
+
&&##$$'$()'$$&$&%%%%%%%(*340++1','',,)-*669:98731/0+&&(%&$$$%(',,4$$,&;>@:2$-')&&25<<>/4.,'*$,--**6<;<>)1...2-&%+///+++*%'*1B:<BFDDB=760998EEA?%33115/5-&&3%%/22:+07660*76..7221680('+058)777(@.A?5.-+&'*7)1*-$)))-/1)/'$*%(**.6*),$&&,+%13-2+'%$(%&+&(&'+%$$+.0023<<6946>93;:8618'/&%,%'()344:=:'*6*5$*39090-/4/4)537-/166+)%,'(%'02)5-5@A:/;75&'3?=9=--=;.,/030'9%&&((48977;<=@;3445%'.2<=;51/(.,),*7*%&(,(()&'/01*)'%%$%%(%%$%$

""")

    # Create adaptor fasta file
    with open(tmp_path / "truseq.fasta", "w") as f:
        f.write(""">truseq_i7
GATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNNNNNATCTCGTATGCCGTCTTCTGCTTG
>truseq_i5
AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT

""")

    fixture = {
        "tmp_path": tmp_path,
        "testdata": tmp_path / "testdata.fastq.gz",
        "index_file": tmp_path / "truseq.fasta",
        "paf": tmp_path / "test.paf",
        "expected_index": "TAACTTGGTC",
    }

    yield fixture

    tmp.cleanup()


def test_run_minimap2(fixture):
    to_test.run_minimap2(
        fastq_in=fixture["testdata"],
        index_file=fixture["index_file"],
        output_paf=fixture["paf"],
        threads=1,
        minimap_b=1,
    )

    expected = "0ad8bdb6-e009-43c5-95b1-d381e699f983\t418\t302\t374\t+\ttruseq_i7\t67\t0\t67\t51\t64\t25\tNM:i:23\tms:i:275\tAS:i:266\tnn:i:10\ttp:A:P\tcm:i:5\ts1:i:29\ts2:i:0\tde:f:0.2388\trl:i:0\tcg:Z:11M2D33M7I21M\tcs:Z::11-ca:6*tg:13*nt*na*na*nc*nt*nt*ng*ng*nt*nc:1*ta:1+ctagaaa:2*gt*tg:17\n0ad8bdb6-e009-43c5-95b1-d381e699f983\t418\t45\t110\t+\ttruseq_i5\t58\t0\t58\t56\t66\t38\tNM:i:10\tms:i:313\tAS:i:305\tnn:i:0\ttp:A:P\tcm:i:10\ts1:i:37\ts2:i:0\tde:f:0.0667\trl:i:0\tcg:Z:15M1D6M7I3M1I33M\tcs:Z::15-a*cg:5+tcccgat:3+g:33\n"
    received = open(fixture["paf"]).read()
    assert expected == received


def test_parse_paf_lines(fixture):
    paf_lines_simple = to_test.parse_paf_lines(fixture["paf"])
    expected_simple = {
        "0ad8bdb6-e009-43c5-95b1-d381e699f983": [
            {
                "read": "0ad8bdb6-e009-43c5-95b1-d381e699f983",
                "adapter": "truseq_i7",
                "rlen": 418,
                "rstart": 302,
                "rend": 374,
                "strand": "+",
                "cg": "cg:Z:11M2D33M7I21M",
                "cs": "cs:Z::11-ca:6*tg:13*nt*na*na*nc*nt*nt*ng*ng*nt*nc:1*ta:1+ctagaaa:2*gt*tg:17",
                "q": 25,
                "iseq": None,
                "sample": None,
            },
            {
                "read": "0ad8bdb6-e009-43c5-95b1-d381e699f983",
                "adapter": "truseq_i5",
                "rlen": 418,
                "rstart": 45,
                "rend": 110,
                "strand": "+",
                "cg": "cg:Z:15M1D6M7I3M1I33M",
                "cs": "cs:Z::15-a*cg:5+tcccgat:3+g:33",
                "q": 38,
                "iseq": None,
                "sample": None,
            },
        ]
    }
    assert paf_lines_simple == expected_simple

    paf_lines_complex = to_test.parse_paf_lines(fixture["paf"])
    expected_complex = {
        "0ad8bdb6-e009-43c5-95b1-d381e699f983": [
            {
                "read": "0ad8bdb6-e009-43c5-95b1-d381e699f983",
                "adapter": "truseq_i7",
                "rlen": 418,
                "rstart": 302,
                "rend": 374,
                "strand": "+",
                "cg": "cg:Z:11M2D33M7I21M",
                "cs": "cs:Z::11-ca:6*tg:13*nt*na*na*nc*nt*nt*ng*ng*nt*nc:1*ta:1+ctagaaa:2*gt*tg:17",
                "q": 25,
                "iseq": None,
                "sample": None,
            },
            {
                "read": "0ad8bdb6-e009-43c5-95b1-d381e699f983",
                "adapter": "truseq_i5",
                "rlen": 418,
                "rstart": 45,
                "rend": 110,
                "strand": "+",
                "cg": "cg:Z:15M1D6M7I3M1I33M",
                "cs": "cs:Z::15-a*cg:5+tcccgat:3+g:33",
                "q": 38,
                "iseq": None,
                "sample": None,
            },
        ]
    }
    assert paf_lines_complex == expected_complex


def test_parse_cd(fixture):
    paf_lines_simple = to_test.parse_paf_lines(fixture["paf"])

    for read_name, alignments in paf_lines_simple.items():
        for alignment in alignments:
            cs_string = alignment["cs"]
            cs_parsed = to_test.parse_cs(
                cs_string=cs_string,
                index=fixture["expected_index"],
                umi_before=0,
                umi_after=0,
            )

            if alignment["adapter"] == "truseq_i7":
                # Perfect match to index
                assert cs_parsed[0] == fixture["expected_index"].lower()
                assert cs_parsed[1] == 0
            elif alignment["adapter"] == "truseq_i5":
                # No index, distance the length of all concatenated substitutions of the cs string
                assert cs_parsed[0] == ""
                assert cs_parsed[1] == 10
            else:
                raise AssertionError("Case not covered.")
