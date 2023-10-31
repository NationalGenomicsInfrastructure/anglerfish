import csv
import Levenshtein as lev
import glob
import re
from dataclasses import dataclass
from itertools import combinations
import yaml
import importlib.resources


p = importlib.resources.files("anglerfish.config").joinpath("adaptors.yaml")
with open(p, "r") as stream:
    adaptors = yaml.safe_load(stream)
delim = "-NNN-"

@dataclass
class SampleSheetEntry:

    sample_name: str
    adaptor: object
    fastq: str
    ont_barcode: str


class Adaptor(object):

    def __init__(self, adaptor, i7_index=None, i5_index=None):

        self.i5 = adaptors[adaptor]["i5"]
        self.i7 = adaptors[adaptor]["i7"]
        self.i5_index = i5_index
        self.i7_index = i7_index
        self.name = f"{adaptor}_len{len(i7_index)}"

        if delim in self.i5 and i5_index is None:
            raise UserWarning("Adaptor has i5 but no sequence was specified")
        if delim in self.i7 and i7_index is None:
            raise UserWarning("Adaptor has i7 but no sequence was specified")

    def get_i5_mask(self):
        if delim in self.i5:
            return self.i5.replace(delim, "N"*len(self.i5_index))
        else:
            return self.i5

    def get_i7_mask(self):
        if delim in self.i7:
            return self.i7.replace(delim, "N"*len(self.i7_index))
        else:
            return self.i7


class SampleSheet(object):

    def __init__(self, input_csv, ont_bc):

        # Read samplesheet in format:
        # sample_name, adaptors, i7_index(-i5_index), fastq_path
        # If we are demuxing a run with ONT barcodes, we will have to assume fastq files are located in "barcode##" folders

        self.samplesheet = []
        try:
            csvfile = open(input_csv, "r")
            dialect = csv.Sniffer().sniff(csvfile.readline(), [',',';','\t'])
            csvfile.seek(0)
            data = csv.DictReader(csvfile,
                fieldnames=['sample_name', 'adaptors', 'index', 'fastq_path'], dialect=dialect)
            rn = 1

            test_globs = {}
            for row in data:
                if row['adaptors'] not in adaptors:
                    raise UserWarning(f"'{row['adaptors']}' not in the list of valid adaptors: {adaptors.keys()}")
                i7i5 = row["index"].split("-")
                i7 = i7i5[0]; i5 = None
                if len(i7i5) > 1:
                    i5 = i7i5[1]

                sample_name = row['sample_name']
                test_globs[row['fastq_path']] = glob.glob(row['fastq_path'])

                bc_re = re.compile("\/(barcode\d\d|unclassified)\/")
                ont_barcode = None
                if ont_bc:
                    ob = re.findall(bc_re, row['fastq_path'])
                    assert len(ob) > 0 and len(ob[0][-1]) > 0, "ONT barcode not found in fastq path. In ONT barcode mode (-n), fastq files must be located in barcode## folders"
                    ont_barcode = ob[0]

                ss_entry = SampleSheetEntry(sample_name, Adaptor(row['adaptors'], i7, i5),row['fastq_path'], ont_barcode)
                self.samplesheet.append(ss_entry)
                rn += 1

            # Explanation: Don't mess around with the globs too much. Don't refer to the same file twice but using globs,
            # e.g, ./input.fastq and ./[i]nput.fastq
            for a,b in combinations(test_globs.values(), 2):
                if len(set(a) & set(b)) > 0:
                    raise UserWarning(f"Fastq paths are inconsistent. Please check samplesheet")

            if not ont_bc and len(set([v[0] for v in test_globs.values()])) > 1:
                raise UserWarning("""Found several different fastq files in samplesheet. Please carefully check any glob patterns. 
                                  If you are using ONT barcodes, please specify the --ont_barcodes flag. Or if you are trying to input several 
                                  sets of fastqs into anglerfish, please run anglerfish separately for each set.""")

        except:
            raise
        finally:
            csvfile.close()


    def minimum_bc_distance(self):
        """ Compute the minimum edit distance between all barcodes in samplesheet, or within each ONT barcode group """

        ss_by_bc = {}
        testset = {}
        for entry in self.samplesheet:
            if entry.ont_barcode in ss_by_bc:
                ss_by_bc[entry.ont_barcode].append(entry.adaptor)
            else:
                ss_by_bc[entry.ont_barcode] = [entry.adaptor]

        for ont_barcode, adaptors in ss_by_bc.items():
            testset[ont_barcode] = []
            for adaptor in adaptors:
                if adaptor.i5_index is not None:
                    testset[ont_barcode].append(adaptor.i5_index+adaptor.i7_index)
                else:
                    testset[ont_barcode].append(adaptor.i7_index)

        fq_distances=[]
        for ont_barcode, adaptors in testset.items():
            distances = []
            if len(adaptors) == 1:
                distances = [len(adaptors[0])]
            else:
                for a, b in [i for i in combinations(adaptors, 2)]:
                    dist = lev.distance(a,b)
                    assert dist > 0, f"""There is one or more identical barcodes in the input samplesheet.
                        First one found: {a}. If these exist in different ONT barcodes, please specify the --ont_barcodes flag."""
                    distances.append(dist)
            fq_distances.append(min(distances))
        return min(fq_distances)

    def get_fastastring(self, adaptor_name=None):

        fastas = {}
        for entry in self.samplesheet:
            if entry.adaptor.name == adaptor_name or adaptor_name is None:
                fastas[entry.adaptor.name+"_i7"] = entry.adaptor.get_i7_mask()
                fastas[entry.adaptor.name+"_i5"] = entry.adaptor.get_i5_mask()

        assert len(fastas) > 0

        outstr = ""
        for key, seq in fastas.items():
            outstr += ">{}\n{}\n".format(key,seq)

        return outstr

    def __iter__(self):
        return iter(self.samplesheet)
    def __next__(self):
        pass
