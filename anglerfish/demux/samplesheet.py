import csv
import glob
import re
from dataclasses import dataclass
from itertools import combinations

import Levenshtein as lev

from anglerfish.demux.adaptor import Adaptor, load_adaptors

ADAPTORS = load_adaptors(raw=True)


@dataclass
class SampleSheetEntry:
    sample_name: str
    adaptor: Adaptor
    fastq: str
    ont_barcode: str | None


class SampleSheet:
    def __init__(self, input_csv: str, ont_barcodes_enabled: bool):
        """Read samplesheet in format:

            sample_name, adaptors, i7_index(-i5_index), fastq_path

        If we are demuxing a run with ONT barcodes, we will have to assume
        fastq files are located in "barcode##" folders.
        """

        self.samplesheet = []
        try:
            csvfile = open(input_csv)
            dialect = csv.Sniffer().sniff(csvfile.readline(), [",", ";", "\t"])
            csvfile.seek(0)
            data = csv.DictReader(
                csvfile,
                fieldnames=["sample_name", "adaptors", "index", "fastq_path"],
                dialect=dialect,
            )
            row_number = 1

            test_globs = {}
            for row in data:
                if row["adaptors"] not in ADAPTORS:
                    raise UserWarning(
                        f"'{row['adaptors']}' not in the list of valid adaptors: {ADAPTORS.keys()}"
                    )
                i7i5_split = row["index"].split("-")
                i7_index = i7i5_split[0]
                if len(i7i5_split) > 1:
                    i5_index = i7i5_split[1]
                else:
                    i5_index = None

                sample_name = row["sample_name"]
                test_globs[row["fastq_path"]] = glob.glob(row["fastq_path"])

                barcode_dir_pattern = re.compile(r"\/(barcode\d\d|unclassified)\/")

                if ont_barcodes_enabled:
                    barcode_dir_match = re.findall(
                        barcode_dir_pattern, row["fastq_path"]
                    )
                    assert (
                        len(barcode_dir_match) > 0 and len(barcode_dir_match[0][-1]) > 0
                    ), "ONT barcode not found in fastq path. In ONT barcode mode (-n), fastq files must be located in barcode## folders"
                    ont_barcode = barcode_dir_match[0]
                else:
                    ont_barcode = None

                ss_entry = SampleSheetEntry(
                    sample_name,
                    Adaptor(
                        name=row["adaptors"],
                        adaptors=ADAPTORS,
                        i5_index=i5_index,
                        i7_index=i7_index,
                    ),
                    row["fastq_path"],
                    ont_barcode,
                )
                self.samplesheet.append(ss_entry)
                row_number += 1

            # Explanation: Don't mess around with the globs too much.
            # Don't refer to the same file twice but using globs, e.g, ./input.fastq and ./[i]nput.fastq
            for a, b in combinations(test_globs.values(), 2):
                if len(set(a) & set(b)) > 0:
                    raise UserWarning(
                        "Fastq paths are inconsistent. Please check samplesheet."
                    )

            if (
                not ont_barcodes_enabled
                and len(set([v[0] for v in test_globs.values()])) > 1
            ):
                raise UserWarning(
                    "Found several different fastq files in samplesheet. Please carefully check any glob patterns."
                    + " If you are using ONT barcodes, please specify the --ont_barcodes flag."
                    + " Or if you are trying to input several sets of fastqs into anglerfish,"
                    + " please run anglerfish separately for each set."
                )

        except:
            raise
        finally:
            csvfile.close()

    def minimum_bc_distance(self) -> int:
        """Compute the minimum edit distance between all barcodes in samplesheet,
        or within each ONT barcode group.
        """

        ont_bc_to_adaptors = {}
        for entry in self.samplesheet:
            if entry.ont_barcode in ont_bc_to_adaptors:
                ont_bc_to_adaptors[entry.ont_barcode].append(entry.adaptor)
            else:
                ont_bc_to_adaptors[entry.ont_barcode] = [entry.adaptor]

        testset = {}
        for ont_barcode, adaptors in ont_bc_to_adaptors.items():
            testset[ont_barcode] = []
            for adaptor in adaptors:
                if adaptor.i5.has_index():
                    testset[ont_barcode].append(adaptor.i5.index + adaptor.i7.index)
                else:
                    testset[ont_barcode].append(adaptor.i7.index)

        min_distances_all_barcodes = []
        for ont_barcode, adaptors in testset.items():
            distances_within_barcode = []
            if len(adaptors) == 1:
                # If there is only one adaptor in the group, the distance is simply the length of the adaptor
                distances_within_barcode = [len(adaptors[0])]
            else:
                for a, b in [i for i in combinations(adaptors, 2)]:
                    distance_this_pair = lev.distance(a, b)
                    assert (
                        distance_this_pair > 0
                    ), f"""There is one or more identical barcodes in the input samplesheet.
                        First one found: {a}. If these exist in different ONT barcodes, please specify the --ont_barcodes flag."""
                    distances_within_barcode.append(distance_this_pair)
            min_distances_all_barcodes.append(min(distances_within_barcode))
        return min(min_distances_all_barcodes)

    def get_fastastring(self, adaptor_name: str = None) -> str:
        fastas = {}
        for entry in self.samplesheet:
            if entry.adaptor.name == adaptor_name or adaptor_name is None:
                fastas[entry.adaptor.name + "_i7"] = entry.adaptor.get_i7_mask()
                fastas[entry.adaptor.name + "_i5"] = entry.adaptor.get_i5_mask()

        assert len(fastas) > 0

        outstr = ""
        for key, seq in fastas.items():
            outstr += f">{key}\n{seq}\n"

        return outstr

    def __iter__(self):
        return iter(self.samplesheet)

    def __next__(self):
        pass
