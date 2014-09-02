import copy
import pprint
from collections import defaultdict
import os.path
import logging

__author__  = 'Soh ISHIGURO'
__email__   = 'si914@sfc.keio.ac.jp'
__license__ = 'Not yet'


class Parser(object):
    
    def __init__(self):
        logging.basicConfig(level=logging.DEBUG,
                            format="%(asctime)s %(levelname)s %(message)s")
    def has_key(self, entry, key):
        return any([attr[0] == key for attr in entry])

    def get_attributes(self, entry, key):
        return filter(lambda attr: attr[0] == key, entry)

    def find_attr(self, entory, query):
        if self.has_key(entory, query):
            return self.get_attributes(entory, query)
        else: return [(query, "")]

    def get_value(self, entry, key):
        return [attr[1] for attr in self.get_attributes(entry, key)]

    def split_line(self, line):
        data = line.split(" - ")
        if len(data) > 1:
            return (data[0], " - ".join(data[1: ]))
        elif line[-2: ] == " -":
            return (line[: -2], "")
        else:
            raise ValueError, "'{}'".format(line)

    def read_ecocyc_file(self, filename):
        if not os.path.isfile(filename):
            raise RuntimeError, "{0} is not found".format(filename)
            
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
                        raise ValueError, "No UNIQUE-ID was specified."
                    retval[unique_id] = copy.copy(cache)
                    unique_id, cache = None, []
                elif line[0] == "/": # Comment line
                    cache[-1] = cache[-1][: -1] + (cache[-1][-1] + line[1: ], )
                elif line[0] == "^": # Flag of modification of avobe attributes
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

    def generate_proteins_entory(self, dat):
        proteins = self.read_ecocyc_file(dat)
        db = defaultdict(dict, {"primary_key": { "source": os.path.basename(dat), "type": "protein"}})
        for protein_id, protein_entry in proteins.items():
            synonyms = self.get_value(protein_entry, "COMMON-NAME")
            if len(synonyms) > 0:
                db.update({protein_id: {"synonyms": synonyms}})

            # COMPONENTS
            components = self.find_attr(protein_entry, "COMPONENTS")
            db[protein_id].update(
                {"COMPONENTS"  : [ (attr[1]) for attr in components],
                 "COEFFICIENTS": [ (1 if len(attr) == 2 else int(attr[3])) for attr in components] })
                
            # CATALYZES
            catalyzes = self.find_attr(protein_entry, "CATALYZES")
            db[protein_id].update(
                {"CATALYZES": [ (_[1]) for _ in catalyzes ] })
             
            # FEATURES
            features = self.find_attr(protein_entry, "FEATURES")
            db[protein_id].update(
                {"FEATURES": [ (_[1]) for _ in features ] })
             
            # GENE
            genes = self.find_attr(protein_entry, "GENE")
            db[protein_id].update(
                {"GENE": [ (_[1]) for _ in genes ] })
                
            # TYPE
            types = self.find_attr(protein_entry, "TYPES")
            db[protein_id].update({"TYPES" : [ (_[1]) for _ in types ] })

            # MODIFIED-FORM
            mod_form = self.find_attr(protein_entry, "MODIFIED-FORM")
            db[protein_id].update({"MODIFIED-FORM" : [ (_[1]) for _ in mod_form ]})

            # UNMODIFIED-FORM
            unmod_form = self.find_attr(protein_entry, "UNMODIFIED-FORM")
            db[protein_id].update({"UNMODIFIED-FORM" : [ (_[1]) for _ in unmod_form ]})

            # REGULATES
            regulates = self.find_attr(protein_entry, "REGULATES")
            db[protein_id].update({"REGULATES" : [ (_[1]) for _ in regulates ]})
            
        return db
        
    def generate_reactions_entory(self, dat):
        reaction = self.read_ecocyc_file(dat)
        db = defaultdict(dict, {"primary_key": {"source": map(os.path.basename, [dat]), "type": "reactions" }})
        
        for reaction_id, reaction_entry in reaction.items():
            # COMMENT
            comment = self.find_attr(reaction_entry, "COMMENT")
            db[reaction_id].update({"COMMENT" : [(_[1]) for _ in comment ]})
            
            # REACTION TYPES                                             
            types = self.find_attr(reaction_entry, "TYPES")              
            db[reaction_id].update({"TYPES" : [(_[1]) for _ in types ]}) 
             
            # ENZYMATIC-REACTION
            enzrec = self.find_attr(reaction_entry, "ENZYMATIC-REACTION")
            db[reaction_id].update({"ENZYMATIC-REACTION": [ (_[1]) for _ in enzrec ]})
             
            # EC-NUMBER
            ec = self.find_attr(reaction_entry, "EC-NUMBER")
            db[reaction_id].update({"EC-NUMBER": [ (_[1]) for _ in ec ]})
             
            # LEFT
            left = self.get_attributes(reaction_entry, "LEFT")
            db[reaction_id].update({"LEFT": [ (_[1]) for _ in left ]})
    
            # RIGHT
            right = self.find_attr(reaction_entry, "RIGHT")
            db[reaction_id].update({"RIGHT": [ (_[1]) for _ in right ]})
             
            # DIRECTION
            direction = self.find_attr(reaction_entry, "REACTION-DIRECTION")
            db[reaction_id].update({"REACTION-DIRECTION": [ (_[1]) for _ in direction ]})
            
        return db

    def generate_features_entory(self, dat):
        return self.read_ecocyc_file(dat) # FIXME
                
    def genrate_enzymes_entory(self, dat):
        return self.read_ecocyc_file(dat) # FIXME


if __name__ == '__main__':
    pass
    
