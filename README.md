## [E-cell split 2014] Reconstruction of protein complex network for whole cell simulation of _E. coli_

### Aim


### Example

```python
# Ecocyc raw dat file
proteins_dat  = '/Users/yukke/dev/ecellp2014/ecocyc/data/proteins.dat'

# generate db
ecoparser = EcocycParser()
proteins_db = ecoparser.generate_proteins_entory(proteins_dat)

# Retrieve all un-modified proteins from protein
mod_p = ModifiedProteinInteraction(proteins_db)
retrieved = mod_p.traceback_to_unmodified_proteins()
....

```
