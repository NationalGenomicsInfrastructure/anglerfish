from __future__ import absolute_import
from __future__ import print_function
import csv
import Levenshtein as lev
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
    }
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

    ## def __iter__(self):
    ## def __next__(self):

    def __init__(self, input_csv):

        # Read samplesheet in format:
        # sample_name, adaptors, i7_index, i5_index

        self.samplesheet = []
        try:
            csvfile = open(input_csv, "r")
            dialect = csv.Sniffer().sniff(csvfile.readline(), [',',';','\t'])
            csvfile.seek(0)
            data = csv.DictReader(csvfile,
                fieldnames=['sample_name', 'adaptors', 'i7_index', 'i5_index'], dialect=dialect)
            for row in data:
                if row['adaptors'] not in adaptors:
                    raise UserWarning("'{}' not in the list of valid adaptors: {}".format(row['adaptors'],adaptors.keys()))
                self.samplesheet.append((row['sample_name'], Adaptor(row['adaptors'], row['i7_index'], row['i5_index'])))
        except:
            raise
        finally:
            csvfile.close()

    def minimum_bc_distance(self):

        testset = []
        for sample, adaptor in self.samplesheet:
            if adaptor.i5_index is not None:
                testset.append(adaptor.i5_index+adaptor.i7_index)
            else:
                testset.append(adaptor.i7_index)

        distances=[]
        for a, b in [i for i in combinations(testset, 2)]:
            distances.append(lev.distance(a,b))

        return min(distances)

    def get_fastastring(self):

        fastas = {}
        for sample, adaptor in self.samplesheet:
            fastas[adaptor.name+"_i7"] = adaptor.get_i7_mask()
            fastas[adaptor.name+"_i5"] = adaptor.get_i5_mask()

        outstr = ""
        for key, seq in fastas.items():
            outstr += ">{}\n{}\n".format(key,seq)

        return outstr
