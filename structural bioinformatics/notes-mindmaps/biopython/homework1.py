# Exercise 1
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
import matplotlib
import matplotlib.pyplot as plt
import copy
import numpy as np
import Bio.PDB
    # Bio.PDB.PDBList
    # Bio.PDB.PDBParser

pdb_id, path = "1ot6", "data/"
Bio.PDB.PDBList().retrieve_pdb_file(pdb_id, pdir=path, file_format='pdb')
structure = Bio.PDB.PDBParser(QUIET = True).get_structure(pdb_id, path + "pdb{}.ent".format(pdb_id))
distances = []
for residue1 in structure[0]['A'] :
    if residue1.id[0] == " " : 
        row = []
        for residue2 in structure[0]['A'] :
            if residue2.id[0] == " " :
                isCB1, isCB2 = False, False
                for atoms in residue1 :
                    if atoms.id == "CB" :   isCB1 = True
                for atoms in residue2 :
                    if atoms.id == "CB" :   isCB2 = True

                if abs(residue1.id[1] - residue2.id[1] >= 1) :
                    if isCB1 and isCB2 :  row.append((int)(residue1["CB"] - residue2["CB"]))
                    if isCB1 and not isCB2 : row.append((int)(residue1["CB"] - residue2["CA"]))
                    if not isCB1 and isCB2 : row.append((int)(residue1["CA"] - residue2["CB"]))
                    if not isCB1 and not isCB2 : row.append((int)(residue1["CA"] - residue2["CA"]))
                else :
                    row.append(None)
    distances.append(row)

finalDistances = np.array(distances, dtype=float)

current_cmap = copy.copy(matplotlib.cm.get_cmap)
#current_cmap.set_bad(color = "white")

fig, ax = plt.subplots(figsize=(12, 12))
im = ax.imshow(finalDistances)
fig.colorbar(im, fraction=0.03, pad=0.05)
plt.savefig(path + "ca_distances_{}.png".format(pdb_id), bbox_inches="tight")

# Plot contact map
contact_map = (finalDistances[:] < 8).astype(float)  
# Calculate the contact map based on a distace threshold 8 Angstrom
fig, ax = plt.subplots(figsize=(12, 12))
im = ax.imshow(contact_map)

# Set ticks
ax.xaxis.set_major_locator(MultipleLocator(10))
ax.xaxis.set_minor_locator(AutoMinorLocator(10))
ax.yaxis.set_major_locator(MultipleLocator(10))
ax.yaxis.set_minor_locator(AutoMinorLocator(10))

plt.savefig(path + 'ca_contacts_{}.png'.format(pdb_id), bbox_inches='tight')


    