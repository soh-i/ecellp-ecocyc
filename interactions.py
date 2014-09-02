from ecocyc_parser import EcocycParser
from operator import add
import itertools
import pprint 

__author__  = 'Soh ISHIGURO'
__email__   = 'si914@sfc.keio.ac.jp'
__license__ = 'Not yet'

pp = pprint.PrettyPrinter()

class ComponentInteraction(object):
    def __init__(self):
        pass

    def parent(self):
        pass
        
class ProteinInteraction(ComponentInteraction):
    def __init__(self):
        ComponentInteraction.__init__(self)

        
class ModifiedProteinInteraction(ComponentInteraction):
    def __init__(self, proteins_db):
        ComponentInteraction.__init__(self)
        self.proteins_db = proteins_db

    def traceback_to_unmodified_proteins(self):
        for p_id in self.proteins_db:
            types = self.proteins_db[p_id].get("TYPES")

            # search "unmodified" protein from "modified" protein
            if types is not None and "Modified-Proteins" in types:
                unmodified_p = self.proteins_db[p_id]["UNMODIFIED-FORM"]
                if unmodified_p[0] in [ pool for pool in self.proteins_db ]:
                    only_type = next(itertools.ifilter(lambda x: "Modified-Proteins" not in x, types), None)
                    print "Modified:'{}' => Original:'{}'".format(only_type, ",".join(unmodified_p))
                    
    def traceback_to_modified_proteins(self):
        pass
                
            
if __name__ == '__main__':
    proteins_dat  = '/Users/yukke/dev/ecellp2014/ecocyc/data/proteins.dat'
    ecoparser = EcocycParser()
    proteins_db = ecoparser.generate_proteins_entory(proteins_dat)
    interaction = ModifiedProteinInteraction(proteins_db)
    print interaction.traceback_to_unmodified_proteins()
