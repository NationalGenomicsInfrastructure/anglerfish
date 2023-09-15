import os
import json
from dataclasses import dataclass
from typing import ClassVar

class Report(object):

    unmatch_header = ["index", "num_reads"]

    def __init__(self, run_name, uuid, version):
        self.run_name = run_name
        self.uuid = uuid
        self.version = version
        self.aln_stats = []
        self.sample_stats = []
        self.unmatched_stats = []

    def add_alignment_stat(self, aln_stat):
        self.aln_stats.append(aln_stat)
    def add_sample_stat(self, sample_stat):
        self.sample_stats.append(sample_stat)
    def add_unmatched_stat(self, unmatched_stat):
        self.unmatched_stats.append(unmatched_stat)

    def write_report(self, outdir):
        with open(os.path.join(outdir,"anglerfish_stats.txt"), "w") as f:
            f.write(f"Anglerfish v. {self.version} (run: {self.run_name}, {self.uuid})\n===================\n")
            for astat in self.aln_stats:
                f.write(f"{astat.adaptor_name}:\n")
                for i,j in astat.paf_stats.items():
                    f.write(f"{j[0]}\t{i} ({j[1]*100:.2f}%)\n")
            f.write("\n{}\n".format("\t".join(getattr(SampleStat, "header"))))
            for sample in self.sample_stats:
                f.write(f"{sample.sample_name}\t{sample.num_reads}\t{sample.mean_read_len}\t{sample.std_read_len}\t{sample.i5_reversed}\n")
            uhead = getattr(Report, 'unmatch_header')
            f.write(f"\n{chr(9).join(uhead)}\n") # chr(9) = tab
            for unmatch in self.unmatched_stats:
                for idx, mnum in unmatch:
                    f.write("{}\t{}\n".format(idx, mnum))
    
    def write_json(self, outdir):
        json_out = {
            "anglerfish_version": self.version,
            "run_name": self.run_name,
            "run_uuid": self.uuid,
            "paf_stats": [],
            "sample_stats": [],
            "undetermined": []
        }
        for astat in self.aln_stats:
            json_out["paf_stats"].append(astat.paf_stats)
        for sample in self.sample_stats:
            slist = [sample.sample_name, sample.num_reads, sample.mean_read_len, sample.std_read_len, sample.i5_reversed]
            json_out["sample_stats"].append(dict(zip(getattr(SampleStat, "header"),slist)))
        for unmatch in self.unmatched_stats:
            for idx, mnum in unmatch:
                json_out["undetermined"].append(dict(zip(getattr(Report, "unmatch_header"),[idx, mnum])))
        with open(os.path.join(outdir,"anglerfish_stats.json"), "w") as f:
            f.write(json.dumps(json_out,indent=2, sort_keys=True))


class AlignmentStat(object):

    def __init__(self, adaptor_name):    
            self.adaptor_name = adaptor_name
            self.paf_stats = {}

    def compute_pafstats(self, num_fq, fragments, singletons, concats, unknowns):
        total = len(fragments)+len(singletons)+len(concats)+len(unknowns)
        self.paf_stats["input_reads"] = [num_fq , 1.0]
        self.paf_stats["reads aligning to adaptor sequences"] = [total, total/float(num_fq)]
        self.paf_stats["aligned reads matching both I7 and I5 adaptor"] = [len(fragments), len(fragments)/float(total)]
        self.paf_stats["aligned reads matching only I7 or I5 adaptor"] = [len(singletons), len(singletons)/float(total)]
        self.paf_stats["aligned reads matching multiple I7/I5 adaptor pairs"] = [len(concats), len(concats)/float(total)]
        self.paf_stats["aligned reads with uncategorized alignments"] = [len(unknowns), len(unknowns)/float(total)]


@dataclass
class SampleStat:

    sample_name: str
    num_reads: int
    mean_read_len: float
    std_read_len: float
    i5_reversed: bool
    ont_barcode: str = None
    header: ClassVar[list] = ["sample_name",
                              "#reads",
                              "mean_read_len",
                              "std_read_len",
                              "i5_reversed"]



