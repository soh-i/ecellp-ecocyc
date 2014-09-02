from operator import add
from collections import defaultdict
import itertools
import pprint

from ecocyc_parser import EcocycParser

__author__  = 'Soh Ishiguro'
__email__   = 'si914@sfc.keio.ac.jp'
__license__ = 'Not yet'

pp = pprint.PrettyPrinter()


class InteractionMap(object):
    def __init__(self):
        pass

    def __exit__(self):
        pass
        
class EnzInteractionMap(InteractionMap):
    def __init__(self):
        InteractionMap.__init__(self)

    def generate_query(self, proteins_db):
        # generate "CATALYZES" from proteins.dat
        queries = list()
        for protein_id in proteins_db:
            cats = proteins_db[protein_id].get("CATALYZES")
            if cats is not None:
                queries.append(cats)
        assert len(cats) < 0, "NO CATALYZES attribtue is found in proteins_db"
            
        if len(queries) > 0:
            return queries
        else: raise ValueError,"{0} is is not contain 'CATALYZES' attribute".format(proteins_db)

    def generate_enz_reaction_map(self, reactions_db="", proteins_db="", debug=False):
        if reactions_db == "": raise ValueError, "Database[{}] name error".format(reactions_db)
        if proteins_db  == "": raise ValueError, "Database[{}] name error".format(proteins_db)
        
        db = defaultdict(dict, {"primary_key": {"source": reactions_db["source"], "type": "enzrxns"}})
        # map(os.path.basename, [proteins_dat, reactions_dat])
        # test_queries = ["UDPNACETYLGLUCOSAMENOLPYRTRANS-ENZRXN", "ENZRXN0-7642", "ENZRXN0-2703"]
        
        queries = self.generate_query(proteins_db)
        for reaction in reactions_db:
            enzrxns = reactions_db[reaction].get("ENZYMATIC-REACTION")
            for query in queries:
                for inn_qry in query:
                    if enzrxns is not None and inn_qry in enzrxns:
                        direction = "".join(map(str, reactions_db[reaction]["REACTION-DIRECTION"]))
                        db.update({inn_qry: {"reaction": reactions_db[reaction], "direction": direction}})
                        if debug:
                            print "ProteinQuery: {}\n EC: {}\n Reactions: {} -> {}, Direction: {}".format(
                                inn_qry, reactions_db[reaction]["EC-NUMBER"],
                                reactions_db[reaction]["LEFT"],
                                reactions_db[reaction]["RIGHT"],
                                reactions_db[reaction]["REACTION-DIRECTION"])
        return db


class ComponentInteraction(object):
    def __init__(self):
        pass

    def parent(self):
        pass
        
        
class ModifiedProteinInteraction(ComponentInteraction):
    def __init__(self, proteins_db):
        ComponentInteraction.__init__(self)
        self.proteins_db = proteins_db
        
    def __repr__(self):
        return ""
        
    def traceback_to_unmodified_proteins(self,):
        mapping = defaultdict(dict, {"primary_key": {"source": ["ecocyc"], "type": "modification" }})
        for p_id in self.proteins_db:
            types = self.proteins_db[p_id].get("TYPES")
            
            # search "unmodified" from "modified"
            if types is not None and "Modified-Proteins" in types:
                unmodified_p = self.proteins_db[p_id]["UNMODIFIED-FORM"]
                assert len(unmodified_p) == 1, "Unexpected elements in unmodified array"
                if unmodified_p[0] in [ _ for _ in self.proteins_db ]:
                    type_name = next(itertools.ifilter(lambda x: "Modified-Proteins" not in x, types), None)
                    mapping.update({unmodified_p[0]: {"TO": type_name, "FROM": unmodified_p[0]}})
        return mapping

                    
if __name__ == '__main__':
    # Raw Ecocyc data
    proteins_dat  = '/Users/yukke/dev/ecellp2014/ecocyc/data/proteins.dat'
    reactions_dat = '/Users/yukke/dev/ecellp2014/ecocyc/data/reactions.dat'
    features_dat  = '/Users/yukke/dev/ecellp2014/ecocyc/data/protein-features.dat'
    enzrxns_dat   = '/Users/yukke/dev/ecellp2014/ecocyc/data/enzrxns.dat'
    
    # Parseing data
    ecoparser    = EcocycParser()
    proteins_db  = ecoparser.generate_proteins_entory(proteins_dat)
    reactions_db = ecoparser.generate_reactions_entory(reactions_dat)
    
    #enz = EnzInteractionMap()
    #enz_interaction = enz.generate_enz_reaction_map(proteins_db=proteins_db, reactions_db=reactions_db)
    #pp.pprint(enz_interaction)

    modproteins = ModifiedProteinInteraction(proteins_db)
    traceback = modproteins.traceback_to_unmodified_proteins()
    for j in traceback:
        print traceback[j]
    
