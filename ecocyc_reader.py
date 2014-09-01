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
        return [attr[1] for attr in self.get_attributes(entry, key)]

    def split_line(self, line):
        data = line.split(" - ")
        if len(data) > 1:
            return (data[0], " - ".join(data[1: ]))
        elif line[-2: ] == " -":
            return (line[: -2], "")
        else:
            raise RuntimeError, "'{}'".format(line)

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
                        raise RuntimeError, "No UNIQUE-ID was specified."
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
        
        db = defaultdict(dict, {"primary_key": { "id": "primary_key", "source": os.path.basename(dat), "type": "protein"}})
        for protein_id, protein_entry in proteins.items():
            synonyms = self.get_value(protein_entry, "COMMON-NAME") # get_value(protein_entry, "SYNONYMS")
            if len(synonyms) > 0:
                db.update({protein_id: {"synonyms": synonyms}})

            # COMPONENTS
            if self.has_key(protein_entry, "COMPONENTS"):
                components = self.get_attributes(protein_entry, "COMPONENTS")
                db["ecocyc"].update(
                    {protein_id: {"COMPONENTS"  : [attr[1] for attr in components],
                                  "COEFFICIENTS": [ (1 if len(attr) == 2 else int(attr[3])) for attr in components] }})
            # CATALYZES
            if self.has_key(protein_entry, "CATALYZES"):
                catalyzes = self.get_attributes(protein_entry, "CATALYZES")
                db["ecocyc"].update(
                    {protein_id: {"CATALYZES": [ (_[1]) for _ in catalyzes ] }})

            # FEATURES
            if self.has_key(protein_entry, "FEATURES"):
                features = self.get_attributes(protein_entry, "FEATURES")
                db["ecocyc"].update(
                    {protein_id: {"FEATURES": [ (_[1]) for _ in features ] }})
            
            # GENE
            if self.has_key(protein_entry, "GENE"):
                genes = self.get_attributes(protein_entry, "GENE")
                db["ecocyc"].update(
                    {protein_id: {"GENE": [ (_[1]) for _ in genes ] }})
        return db

    def generate_features_entory(self, dat):
        features = self.read_ecocyc_file(dat)
        return features
        
    def generate_reactions_entory(self, dat):
        reaction = self.read_ecocyc_file(dat)
        db = defaultdict(dict, {"primary_key":
                                { "id": "primary_key", "source": os.path.basename(dat), "type": "reaction" }})
        
        def find_attr(entory, query):
            if self.has_key(entory, query):
                return self.get_attributes(entory, query)
            else: return [None]
            
        for reaction_id, reaction_entry in reaction.items():
            # REACTION TYPES
            if self.has_key(reaction_entry, "TYPES"):
                types = self.get_attributes(reaction_entry, "TYPES")
                db["ecocyc"].update( {reaction_id: {"TYPES" : [(_[1]) for _ in types ] }})

            # ENZYMATIC-REACTION
            if self.has_key(reaction_entry, "ENZYMATIC-REACTION"):
                enzrec = self.get_attributes(reaction_entry, "ENZYMATIC-REACTION")
                db["ENZYMATIC-REACTION"].update(
                    {reaction_id: {"ENZYMATIC-REACTION": [ (_[1]) for _ in enzrec ] }})

            # EC-NUMBER
            if self.has_key(reaction_entry, "EC-NUMBER"):
                ec = self.get_attributes(reaction_entry, "EC-NUMBER")
                db["ecocyc"].update( {reaction_id: {"EC-NUMBER": [ (_[1]) for _ in ec ] }})

            # LEFT
            left = find_attr(reaction_entry, "LEFT")
            db["ecocyc"].update( {reaction_id: {"LEFT": [ (_[1]) for _ in left ] }})
                
            # RIGHT
            right = find_attr(reaction_entry, "RIGHT")
            db["ecocyc"].update( {reaction_id: {"RIGHT": [ (_[1]) for _ in right ] }})

            # DIRECTION
            direction = find_attr(reaction_entry, "REACTION-DIRECTION")
            db["ecocyc"].update( {reaction_id: {"REACTION-DIRECTION": [ (_[1]) for _ in direction ] }})
                        
        return db

    def genrate_enzymes_entory(self, dat):
        enzyme = self.read_ecocyc_file(dat)
        return enzyme

        
class EcocycQuery(object):
    def __init__(self):
        pass

    def query():
        pass
        
        
if __name__ == '__main__':
    proteins_dat = '/Users/yukke/dev/ecellp2014/ecocyc/data/proteins.dat'
    features_dat = '/Users/yukke/dev/ecellp2014/ecocyc/data/protein-features.dat'
    reactions_dat = '/Users/yukke/dev/ecellp2014/ecocyc/data/reactions.dat'
    enzrxns_dat = '/Users/yukke/dev/ecellp2014/ecocyc/data/enzrxns.dat'
    
    ecoparser = EcocycParser()
    #proteins_db = ecoparser.generate_proteins_entory(proteins_dat)
    #features_db = ecoparser.generate_features_entory(features_dat)
    reac = ecoparser.generate_reactions_entory(reactions_dat)
    print reac['ecocyc']

    
    
    #print features_db.keys()[:10]
    #pp.pprint(protein_db)
    
