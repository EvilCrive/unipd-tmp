import Bio.PDB
import Bio.PDB.PDBParser
import Bio.SeqUtils
import Bio.PDB.PDBIO
import Bio.SeqIO.PdbIO

#fetch pdb file 1ot6
pdb_id = '1ot6'
path = 'data/'
pdbL = Bio.PDB.PDBList()
pdbL.retrieve_pdb_file(pdb_id, pdir = path, file_format = 'pdb')

# get seqres, iterating for each chain (in this case only 1)
with open(path + "pdb{}.ent".format(pdb_id)) as f :
    seq_records = (Bio.SeqIO.PdbIO.PdbSeqresIterator(f))
    for s in seq_records :
        print(s)

#load structure
structure = Bio.PDB.PDBParser(QUIET = True)
structure = structure.get_structure(pdb_id, path+"pdb{}.ent".format(pdb_id))

for model in structure :
    count = 0
    for chain in model :
        for residue in chain :
            if not Bio.PDB.is_aa(residue) :
                print("model {} chain {} residue_id {} resname {} resname_3to1 {}".format(model.id, chain.id, residue.id, residue.get_resname(), Bio.SeqUtils.IUPACData.protein_letters_3to1.get(residue.get_resname().capitalize())))

# extract list of residues
domain_residues = []
start_flag = False
domain_start = (" ", 10, " ")
domain_end = (" ", 100, " ")
for residue in structure[0]['A'].get_residues() :
    if residue.id[0] == " " :
        if residue.id == domain_start :
            start_flag = True
        if start_flag :
            domain_residues.append(residue)
        if residue.id == domain_end :
            break
print(domain_residues)

#save a pdb chain
class Select(Bio.PDB.Select) :
    def __init__(self, chain_ids = None, residues = None):
        self.chain_ids = chain_ids
        self.residues = residues
    def accept_chain(self, chain) :
        return self.chain_ids is None or chain.id in self.chain_ids

    def accept_residue(self, residue):
        return self.residues is None or residue in self.residues

    def accept_atom(self, atom):
        return not atom.is_disordered() or atom.get_altloc() == "A"