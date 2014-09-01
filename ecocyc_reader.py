import copy
import pprint
from collections import defaultdict
import os.path

pp = pprint.PrettyPrinter()


class Parser(object):
    def __init__(self):
        pass

    def has_key(self, entry, key):
        return any([attr[0] == key for attr in entry])

    def get_attributes(self, entry, key):
        return filter(lambda attr: attr[0] == key, entry)

    def get_value(self, entry, key):
        return [attr[1] for attr in get_attributes(entry, key)]

    def split_line(self, line):
        data = line.split(" - ")
        if len(data) > 1:
            return (data[0], " - ".join(data[1: ]))
        elif line[-2: ] == " -":
            return (line[: -2], "")
        else:
            raise RuntimeError, "'{}'".format(line)

    def read_ecocyc_file(self, filename):
        unique_id, cache = None, []
        retval = {}
        with open(filename, "r") as fin:
            for line in fin:
                if line is None or line == "":
                    break
                line = line.strip()
                if line[0] == "#":
                    pass # Header line
                elif line == "//":
                    if unique_id is None:
                        raise RuntimeError, "No UNIQUE-ID was specified."
                    retval[unique_id] = copy.copy(cache)
                    unique_id, cache = None, []
                elif line[0] == "/":
                    cache[-1] = cache[-1][: -1] + (cache[-1][-1] + line[1: ], )
                elif line[0] == "^":
                    cache[-1] = cache[-1] + self.split_line(line[1: ])
                else:
                    key, value = self.split_line(line)
                    if key == "UNIQUE-ID":
                        unique_id = value
                    cache.append(self.split_line(line))
        return retval

        
class EcocycParser(Parser):
    def __init__(self):
        Parser.__init__(self)
                
    def read_ecocyc_file(self, filename):
        unique_id, cache = None, []
        retval = {}
        with open(filename, "r") as fin:
            for line in fin:
                if line is None or line == "":
                    break
                line = line.strip()
                if line[0] == "#":
                    pass # Header line
                elif line == "//":
                    if unique_id is None:
                        raise RuntimeError, "No UNIQUE-ID was specified."
                    retval[unique_id] = copy.copy(cache)
                    unique_id, cache = None, []
                elif line[0] == "/":
                    cache[-1] = cache[-1][: -1] + (cache[-1][-1] + line[1: ], )
                elif line[0] == "^":
                    cache[-1] = cache[-1] + self.split_line(line[1: ])
                else:
                    key, value = self.split_line(line)
                    if key == "UNIQUE-ID":
                        unique_id = value
                    cache.append(self.split_line(line))
        return retval

    def generate_protein_entory(self):
        for protein_id, protein_entry in proteins.items():
            doc = defaultdict(dict, {"id": protein_id, "source": os.path.basename(proteins_data), "type": "protein"})
            synonyms = get_value(protein_entry, "COMMON-NAME") # get_value(protein_entry, "SYNONYMS")
        
            if len(synonyms) > 0:
                doc.update({"synonyms": synonyms})
            
            if has_key(protein_entry, "COMPONENTS"):
                components = get_attributes(protein_entry, "COMPONENTS")
                doc["ecocyc"].update(
                    {"COMPONENTS": [attr[1] for attr in components],
                     "COEFFICIENTS": [ (1 if len(attr) == 2 else int(attr[3])) for attr in components] })
            
            if has_key(protein_entry, "CATALYZES"):
                catalyzes = get_attributes(protein_entry, "CATALYZES")
                doc["ecocyc"].update({"CATALYZES": [ (_[1]) for _ in catalyzes ] })

            if has_key(protein_entry, "FEATURES"):
                features = get_attributes(protein_entry, "FEATURES")
                doc["ecocyc"].update({"FEATURES": [ (_[1]) for _ in features ] })

            if has_key(protein_entry, "GENE"):
                genes = get_attributes(protein_entry, "GENE")
                doc["ecocyc"].update({"GENE": [ (_[1]) for _ in genes ] })
        return doc


if __name__ == '__main__':
    proteins_data = '/Users/yukke/dev/ecellp2014/ecocyc/data/proteins.dat'
    parser = EcocycParser()
    protein_db = parser.read_ecocyc_file(proteins_data)
    print protein_db['EG12874-MONOMER']
    
