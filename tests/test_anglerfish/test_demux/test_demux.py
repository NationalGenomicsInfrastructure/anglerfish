import tempfile
from pathlib import Path

from anglerfish.demux import demux as to_test


def test_run_minimap2():
    tmp = tempfile.TemporaryDirectory()
    tmp_path = Path(tmp.name)
    paf_path = tmp_path / "test.paf"

    to_test.run_minimap2(
        fastq_in="testdata/explore/explore_testdata.fastq",
        index_file="testdata/explore/illumina_ud.fasta",
        output_paf=paf_path,
        threads=1,
        minimap_b=1,
    )

    expected = (
        "\n".join(
            [
                "truseq_i5	58	0	58	+	truseq-read_10_C-index_5G	130	0	58	58	58	60	NM:i:0	ms:i:348	AS:i:348	nn:i:0	tp:A:P	cm:i:16	s1:i:54	s2:i:0	de:f:0	rl:i:0	cg:Z:58M	cs:Z::58",
                "truseq_i7	57	0	57	+	truseq-read_10_C-index_5G	130	68	130	57	62	60	NM:i:5	ms:i:333	AS:i:328	nn:i:0	tp:A:P	cm:i:15	s1:i:54	s2:i:0	de:f:0.0172	rl:i:0	cg:Z:33M5D24M	cs:Z::33-ggggg:24",
            ]
        )
        + "\n"
    )

    assert expected == open(paf_path).read()