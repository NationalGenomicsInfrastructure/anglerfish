from __future__ import absolute_import
from __future__ import print_function
import csv
import Levenshtein as lev
import os
from itertools import combinations

adaptors = {
    "truseq": {
                "i5": ["AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT"],
                "i7": ["GATCGGAAGAGCACACGTCTGAACTCCAGTCAC",True,"ATCTCGTATGCCGTCTTCTGCTTG"]
    },
    "truseq_dual": {
                "i5": ["AATGATACGGCGACCACCGAGATCTACAC",True,"ACACTCTTTCCCTACACGACGCTCTTCCGATCT"],
                "i7": ["GATCGGAAGAGCACACGTCTGAACTCCAGTCAC",True,"ATCTCGTATGCCGTCTTCTGCTTG"]
    },
    "nextera_legacy": {
                "i5": ["AATGATACGGCGACCACCGAGATCTACACGCCTCCCTCGCGCCATCAG"],
                "i7": ["CAAGCAGAAGACGGCATACGAGAT",True,"CGGTCTGCCTTGCCAGCCCGCTCAG"]
    },
    "nextera_dual": {
                "i5": ["AATGATACGGCGACCACCGAGATCTACAC",True,"GTCTCGTGGGCTCGG"],
                "i7": ["CAAGCAGAAGACGGCATACGAGAT",True,"ATCTCGTATGCCGTCTTCTGCTTG"]
    },
}
class Adaptor(object):

    def __init__(self, adaptor, i7_index=None, i5_index=None):

        self.i5_list = adaptors[adaptor]["i5"]
        self.i7_list = adaptors[adaptor]["i7"]
        self.i5_index = i5_index
        self.i7_index = i7_index
        self.name = "{}_len{}".format(adaptor, len(i7_index))

        if True in self.i5_list and i5_index is None:
            raise UserWarning("Adaptor has i5 but no sequence was specified")
        if True in self.i7_list and i7_index is None:
            raise UserWarning("Adaptor has i7 but no sequence was specified")

    def get_i5_mask(self):
        if True in self.i5_list:
            return "".join(map(lambda x: "".join(["N"]*len(self.i5_index)) if x==True else x, self.i5_list))
        else:
            return "".join(self.i5_list)

    def get_i7_mask(self):
        if True in self.i7_list:
            return "".join(map(lambda x: "".join(["N"]*len(self.i7_index)) if x==True else x, self.i7_list))
        else:
            return "".join(self.i7_list)



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
            for row in data:
                if row['adaptors'] not in adaptors:
                    raise UserWarning("'{}' not in the list of valid adaptors: {}".format(row['adaptors'],adaptors.keys()))
                i7i5 = row["index"].split("-")
                i7 = i7i5[0]; i5 = None
                if len(i7i5) > 1:
                    i5 = i7i5[1]

                sample_name = row['sample_name']
                # TODO: find a more clever way of resolving duplicate names
                if row['sample_name'] in [sn[0] for sn in self.samplesheet]:
                    sample_name = sample_name+"_row"+str(rn)
                assert os.path.exists(row['fastq_path']) == 1
                self.samplesheet.append((sample_name, Adaptor(row['adaptors'], i7, i5),row['fastq_path']))
                rn += 1
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
