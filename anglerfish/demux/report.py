import json
import logging
import os
from dataclasses import asdict, dataclass
from typing import ClassVar

log = logging.getLogger("anglerfish")


class Report:
    unmatch_header = ["index", "num_reads", "ont_barcode"]

    def __init__(self, run_name, uuid, version):
        self.run_name = run_name
        self.uuid = uuid
        self.version = version
        self.aln_stats = []
        self.sample_stats = []
        self.unmatched_stats = {}

    def add_alignment_stat(self, aln_stat):
        self.aln_stats.append(aln_stat)

    def add_sample_stat(self, sample_stat):
        self.sample_stats.append(sample_stat)

    def add_unmatched_stat(self, unmatched_stat, ont_barcode, adaptor_name):
        self.unmatched_stats[(ont_barcode, adaptor_name)] = unmatched_stat

    def write_report(self, outdir):
        with open(os.path.join(outdir, "anglerfish_stats.txt"), "w") as f:
            f.write(
                f"Anglerfish v. {self.version} (run: {self.run_name}, {self.uuid})\n===================\n"
            )
            for astat in self.aln_stats:
                f.write(f"{astat.adaptor_name}:\n")
                for i, j in astat.paf_stats.items():
                    f.write(f"{j[0]}\t{i} ({j[1]*100:.2f}%)\n")
            f.write("\n{}\n".format("\t".join(getattr(SampleStat, "header"))))
            for sample in self.sample_stats:
                f.write(
                    f"{sample.sample_name}\t{sample.num_reads}\t{sample.mean_read_len}\t{sample.std_read_len}\t{sample.i7_reversed}\t{sample.i5_reversed}\t{sample.ont_barcode}\n"
                )
            uhead = getattr(Report, "unmatch_header")
            f.write(f"\n{chr(9).join(uhead)}\n")  # chr(9) = tab
            for key, unmatch in self.unmatched_stats.items():
                for idx, mnum in unmatch:
                    f.write(f"{idx}\t{mnum}\t{key[0]}\n")
        log.debug(
            f"Wrote anglerfish_stats.txt to {outdir}, size: {os.path.getsize(os.path.join(outdir, 'anglerfish_stats.txt'))} bytes"
        )

    def write_json(self, outdir):
        json_out = {
            "anglerfish_version": self.version,
            "run_name": self.run_name,
            "run_uuid": self.uuid,
            "paf_stats": [],
            "sample_stats": [],
            "undetermined": [],
        }
        for astat in self.aln_stats:
            json_out["paf_stats"].append(astat.paf_stats)
        for sample in self.sample_stats:
            slist = [
                sample.sample_name,
                sample.num_reads,
                sample.mean_read_len,
                sample.std_read_len,
                sample.i7_reversed,
                sample.i5_reversed,
                sample.ont_barcode,
            ]
            json_out["sample_stats"].append(
                dict(zip(getattr(SampleStat, "header"), slist))
            )
        for key, unmatch in self.unmatched_stats.items():
            for idx, mnum in unmatch:
                json_out["undetermined"].append(
                    dict(zip(getattr(Report, "unmatch_header"), [idx, mnum, key[0]]))
                )
        with open(os.path.join(outdir, "anglerfish_stats.json"), "w") as f:
            f.write(json.dumps(json_out, indent=2, sort_keys=True))
            log.debug(
                f"Wrote anglerfish_stats.json to {outdir}, size: {os.path.getsize(os.path.join(outdir, 'anglerfish_stats.json'))} bytes"
            )

    def write_dataframe(self, outdir, samplesheet):
        """Write a dataframe of the stats to a csv file.
        TODO: This needs be cleaned up and made more robust. Especially lock in / decouple from upstream the header names and order:
        sample_name, num_reads, mean_read_len, std_read_len, i7_reversed, i5_reversed, ont_barcode, adaptor_name, i7_index, i5_index
        """
        out_list = []
        for sample in self.sample_stats:
            s_dict = asdict(sample)
            for sentry in samplesheet:
                sen_dict = asdict(sentry)
                if sen_dict["sample_name"] == s_dict["sample_name"]:
                    s_dict["adaptor_name"] = sen_dict["adaptor"].name
                    s_dict["i7_index"] = sen_dict["adaptor"].i7.index
                    s_dict["i5_index"] = sen_dict["adaptor"].i5.index
            out_list.append(s_dict)
        for key, unmatch in self.unmatched_stats.items():
            for unmatch_sample in unmatch:
                un = {i: None for i in out_list[-1].keys()}
                i7i5 = [i.upper() for i in unmatch_sample[0].split("+")]
                if len(i7i5) == 1:
                    i7i5.append(None)
                un["adaptor_name"] = key[1]
                un["num_reads"] = unmatch_sample[1]
                un["ont_barcode"] = key[0]
                un["i7_index"] = i7i5[0]
                un["i5_index"] = i7i5[1]
                out_list.append(un)
        with open(os.path.join(outdir, "anglerfish_dataframe.csv"), "w") as f:
            out_header = out_list[0].keys()
            f.write(",".join(out_header))
            f.write("\n")
            for out in out_list:
                f.write(",".join([str(out[i]) for i in out_header]))
                f.write("\n")
        log.debug(
            f"Wrote anglerfish_dataframe.csv to {outdir}, size: {os.path.getsize(os.path.join(outdir, 'anglerfish_dataframe.csv'))} bytes"
        )


class AlignmentStat:
    def __init__(self, adaptor_name):
        self.adaptor_name = adaptor_name
        self.paf_stats = {}

    def compute_pafstats(self, num_fq, fragments, singletons, concats, unknowns):
        total = len(fragments) + len(singletons) + len(concats) + len(unknowns)
        self.paf_stats["input_reads"] = [num_fq, 1.0]
        self.paf_stats["reads aligning to adaptor sequences"] = [
            total,
            total / float(num_fq),
        ]
        self.paf_stats["aligned reads matching both I7 and I5 adaptor"] = [
            len(fragments),
            len(fragments) / float(total),
        ]
        self.paf_stats["aligned reads matching only I7 or I5 adaptor"] = [
            len(singletons),
            len(singletons) / float(total),
        ]
        self.paf_stats["aligned reads matching multiple I7/I5 adaptor pairs"] = [
            len(concats),
            len(concats) / float(total),
        ]
        self.paf_stats["aligned reads with uncategorized alignments"] = [
            len(unknowns),
            len(unknowns) / float(total),
        ]


@dataclass
class SampleStat:
    sample_name: str
    num_reads: int
    mean_read_len: float
    std_read_len: float
    i7_reversed: bool
    i5_reversed: bool
    ont_barcode: str | None = None
    header: ClassVar[list] = [
        "sample_name",
        "#reads",  # We specify this for historical reasons
        "mean_read_len",
        "std_read_len",
        "i7_reversed",
        "i5_reversed",
        "ont_barcode",
    ]
