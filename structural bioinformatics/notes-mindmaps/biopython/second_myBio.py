import copy
from locale import currency
import pdb
from Bio.PDB import PDBList, NeighborSearch
from Bio.PDB.PDBParser import PDBParser
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
import scipy.signal

def get_distance_matrix(residues, seq_sep = 6) :
    distances = []
    for residue1 in residues :
        if residue1.id[0] == " " :
            row = []
            for residue2 in residues :
                if residue2.id[0] == " " :
                    if abs(residue1.id[1] - residue2.id[1] >= seq_sep) :
                        row.append(residue1["CA"] - residue2["CA"])
                    else :
                        row.append(None)
            distances.append(row)
    return np.array(distances, dtype=float)

pdb_id, path = "1ucd", "data/"
PDBList().retrieve_pdb_file(pdb_id, pdir=path, file_format='pdb')
structure = PDBParser(QUIET = True).get_structure(pdb_id, path + "pdb{}.ent".format(pdb_id))

dist_matrix = get_distance_matrix(structure[0]['A'], 6)

current_cmap = copy.copy(matplotlib.cm.get:)    


