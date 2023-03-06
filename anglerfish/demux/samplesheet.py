from __future__ import absolute_import
from __future__ import print_function
import csv
import Levenshtein as lev
import glob
from itertools import combinations
import yaml
import importlib.resources


p = importlib.resources.files("anglerfish.config").joinpath("adaptors.yaml")
with open(p, "r") as stream:
    adaptors = yaml.safe_load(stream)
delim = "-NNN-"


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

    def __init__(self, input_csv):

        # Read samplesheet in format:
        # sample_name, adaptors, i7_index(-i5_index), fastq_path

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
                if row['sample_name'] in [sn[0] for sn in self.samplesheet]:
                    sample_name = sample_name+"_row"+str(rn)
                assert len(glob.glob(row['fastq_path'])) > 0
                test_globs[row['fastq_path']] = glob.glob(row['fastq_path'])
                self.samplesheet.append((sample_name, Adaptor(row['adaptors'], i7, i5),row['fastq_path']))
                rn += 1

            # Explaination: Don't mess around with the globs too much.
            for a,b in combinations(test_globs.values(), 2):
                if len(set(a) & set(b)) > 0:
                    raise UserWarning(f"Fastq paths are inconsistent. Please check samplesheet")
        except:
            raise
        finally:
            csvfile.close()

    def minimum_bc_distance(self):

        ss_by_fastq = {}
        testset = {}
        for _, adaptor, fastq in self.samplesheet:
            if fastq in ss_by_fastq:
                ss_by_fastq[fastq].append(adaptor)
            else:
                ss_by_fastq[fastq] = [adaptor]

        for fastq, adaptors in ss_by_fastq.items():
            testset[fastq] = []
            for adaptor in adaptors:
                if adaptor.i5_index is not None:
                    testset[fastq].append(adaptor.i5_index+adaptor.i7_index)
                else:
                    testset[fastq].append(adaptor.i7_index)

        fq_distances=[]
        for fastq, adaptors in testset.items():
            distances = []
            if len(adaptors) == 1:
                distances = [len(adaptors[0])]
            else:
                for a, b in [i for i in combinations(adaptors, 2)]:
                    distances.append(lev.distance(a,b))
            fq_distances.append(min(distances))
        return min(fq_distances)

    def get_fastastring(self, adaptor_name=None):

        fastas = {}
        for sample, adaptor, fastq in self.samplesheet:
            if adaptor.name == adaptor_name or adaptor_name is None:
                fastas[adaptor.name+"_i7"] = adaptor.get_i7_mask()
                fastas[adaptor.name+"_i5"] = adaptor.get_i5_mask()

        assert len(fastas) > 0

        outstr = ""
        for key, seq in fastas.items():
            outstr += ">{}\n{}\n".format(key,seq)

        return outstr

    def __iter__(self):
        return iter(self.samplesheet)
    def __next__(self):
        pass
