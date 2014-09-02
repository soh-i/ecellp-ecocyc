import copy
import pprint
from collections import defaultdict
import os.path
import logging

__author__  = 'Soh ISHIGURO'
__email__   = 'si914@sfc.keio.ac.jp'
__license__ = 'Not yet'

pp = pprint.PrettyPrinter()
logging.basicConfig(level=logging.DEBUG, format="%(asctime)s %(levelname)s %(message)s")


class Parser(object):
    def __init__(self):
        pass

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


        
class InteractionMap(object):
    def __init__(self):
        pass

    def __exit__(self):
        pass

class EnzInteractionMap(InteractionMap):
    def __init__(self):
        InteractionMap.__init__(self)

    def generate_query(self, proteins_db):
        # search catalyzes attribute in proteins.dat
        queries = list()
        for protein_id in proteins_db:
            cats = proteins_db[protein_id].get("CATALYZES")
            if cats is not None:
                queries.append(cats)
        if len(queries) > 0:
            return queries
        else: raise ValueError,"{0} is is not contain 'CATALYZES' attribute".format(proteins_db)

    def generate_enz_reaction_map(self, reactions_db="", proteins_db="", debug=False):
        if reactions_db == "": raise ValueError, "Database[{}] name error".format(reactions_db)
        if proteins_db  == "": raise ValueError, "Database[{}] name error".format(proteins_db)
        
        db = defaultdict(dict, {"primary_key": {"source": reactions_db["source"], "type": "enzrxns"}})
        # map(os.path.basename, [proteins_dat, reactions_dat])
        # proteins_query_t = ["UDPNACETYLGLUCOSAMENOLPYRTRANS-ENZRXN", "ENZRXN0-7642", "ENZRXN0-2703"]
        
        queries = self.generate_query(proteins_db)
        for reaction in reactions_db:
            enzrxns = reactions_db[reaction].get("ENZYMATIC-REACTION")
            for query in queries[:]:
                #for query in [ proteins_db[protein_id].get("CATALYZES")
                #for protein_id in proteins_db if proteins_db[protein_id].get("CATALYZES") is not None]:
                for inn_qry in query:
                    if enzrxns is not None and inn_qry in enzrxns:
                        direction = "".join(map(str, reactions_db[reaction]["REACTION-DIRECTION"]))
                        db.update({inn_qry: {"reaction": reactions_db[reaction], "direction": direction}})
                        
                        if debug:
                            logging.debug("ProteinQuery: {}\n EC: {}\n Reactions: {} -> {}, Direction: {}".format(
                                inn_qry, reactions_db[reaction]["EC-NUMBER"],
                                reactions_db[reaction]["LEFT"],
                                reactions_db[reaction]["RIGHT"],
                                reactions_db[reaction]["REACTION-DIRECTION"]))
        return db


if __name__ == '__main__':
    proteins_dat  = '/Users/yukke/dev/ecellp2014/ecocyc/data/proteins.dat'
    features_dat  = '/Users/yukke/dev/ecellp2014/ecocyc/data/protein-features.dat'
    reactions_dat = '/Users/yukke/dev/ecellp2014/ecocyc/data/reactions.dat'
    enzrxns_dat   = '/Users/yukke/dev/ecellp2014/ecocyc/data/enzrxns.dat'
    
    ecoparser = EcocycParser()
    proteins_db = ecoparser.generate_proteins_entory(proteins_dat)
    reactions_db = ecoparser.generate_reactions_entory(reactions_dat)

    #rec_map = EnzInteractionMap()
    #all_enzyme_reactions = rec_map.generate_enz_reaction_map(proteins_db=proteins_db,
    #                                                         reactions_db=reactions_db)

    
    # query: proteins, db: reactions
    #for rec in reactions_db:
    #    enzrecs = reactions_db[rec].get("ENZYMATIC-REACTION")
    #    if enzrecs is not None:
    #        pass
    #        
    #        #if "THYKI-ENZRXN" in enzrec for r in recs:
    #        #print reactions_db[rec]
    #        #if "ENZRXN0-1621" in enzrec:
    #        #print  reactions_db[rec].get("RIGHT")
    #        #print reactions_db[rec].get("LEFT")
    #        #print reactions_db[rec].get("REACTION-DIRECTION")

    
