from anglerfish.demux import demux as to_test


def test_run_minimap2(demux_fixture):
    """Check that the function runs successfully, not the output."""
    # Test alignment on single read
    to_test.run_minimap2(
        fastq_in=demux_fixture["testdata_single"],
        index_file=demux_fixture["index_file"],
        output_paf=demux_fixture["paf_single"],
        threads=1,
        minimap_b=1,
    )

    # Create aligntment from multiple reads
    to_test.run_minimap2(
        fastq_in=demux_fixture["testdata_multiple"],
        index_file=demux_fixture["index_file"],
        output_paf=demux_fixture["paf_multiple"],
        threads=1,
        minimap_b=1,
    )


def test_map_reads_to_alns(demux_fixture):
    reads_alns = to_test.map_reads_to_alns(demux_fixture["paf_single"])

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


def test_parse_cs(demux_fixture):
    test_cs_str = (
        "cs:Z::11-ca:6*tg:13*nt*na*na*nc*nt*nt*ng*ng*nt*nc:1*ta:1+ctagaaa:2*gt*tg:17"
    )
    expected_cs_parsed = ("taacttggtc", 0)

    cs_parsed = to_test.parse_cs(
        cs_string=test_cs_str,
        index_seq=demux_fixture["index"],
        umi_before=0,
        umi_after=0,
    )

    assert cs_parsed == expected_cs_parsed


def test_categorize_matches(demux_fixture):
    i5_name = "truseq_i5"
    i7_name = "truseq_i7"
    reads_alns = to_test.map_reads_to_alns(demux_fixture["paf_multiple"])

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
