from ecocyc_parser import EcocycParser

__author__  = 'Soh ISHIGURO'
__email__   = 'si914@sfc.keio.ac.jp'
__license__ = 'Not yet'


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

    def modified_proteins(self):
        for p_id in self.proteins_db:
            #self.proteins_db[p_id].get("UNMODIFIED-FORM")
            #print self.proteins_db[p_id]
            types =  self.proteins_db[p_id].get("TYPES") # == "Modified-Proteins"
            print types
            #print self.proteins_db[p_id]
            #print p_id

            
if __name__ == '__main__':
    proteins_dat  = '/Users/yukke/dev/ecellp2014/ecocyc/data/proteins.dat'
    ecoparser = EcocycParser()
    proteins_db = ecoparser.generate_proteins_entory(proteins_dat)
    interaction = ModifiedProteinInteraction(proteins_db)
    print interaction.modified_proteins()


    
