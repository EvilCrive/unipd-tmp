import re
from Bio.PDB import PDBList, is_aa, PDBIO, NeighborSearch
from Bio.PDB.PDBParser import PDBParser
from Bio.SeqUtils import IUPACData, seq1
from Bio.PDB.PDBIO import Select
from Bio.SeqIO.PdbIO import PdbSeqresIterator
import copy
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
import scipy.signal

pdb_id = "2k5d"
chain_id = "A"
path = "data\\"

m_matrix = [
    [-0.20, -0.44, 0.16, 0.26, -0.46, -0.26, 0.50, -0.57, 0.10, -0.36, -0.22, 0.07, 0.14, 0.01, 0.20, -0.09, -0.05, -0.42, 0.05, -0.50],
    [-0.44, -2.99, 0.21, 0.19, -0.88, -0.34, -1.11, -0.36, 0.09, -0.53, -0.43, -0.52, -0.14, -0.43, -0.24, 0.13, -0.22, -0.62, 0.24, -0.79],
    [0.16, 0.21, 0.17, 0.55, 0.38, 0.35, -0.23, 0.44, 0.39, 0.28, 0.35, -0.02, 1.03, 0.49, -0.37, 0.19, -0.12, 0.69, 0.04, 0.43],
    [0.26, 0.19, 0.55, 0.60, 0.55, 0.65, 0.18, 0.37, 0.47, 0.33, 0.29, 0.01, 0.69, 0.04, -0.52, 0.18, 0.37, 0.39, 0.03, 0.17],
    [-0.46, -0.88, 0.38, 0.55, -0.94, 0.17, -0.40, -0.88, 0.01, -1.08, -0.78, 0.22, 0.20, 0.26, -0.19, -0.22, 0.02, -1.15, -0.60, -0.88],
    [-0.26, -0.34, 0.35, 0.65, 0.17, -0.12, 0.18, 0.24, 0.19, 0.24, 0.02, -0.04, 0.60, 0.46, 0.50, 0.28, 0.28, 0.27, 0.51, -0.35],
    [0.50, -1.11, -0.23, 0.18, -0.40, 0.18, 0.42, -0.00, 0.79, -0.24, -0.07, 0.20, 0.25, 0.69, 0.24, 0.21, 0.11, 0.16, -0.85, -0.26],
    [-0.57, -0.36, 0.44, 0.37, -0.88, 0.24, -0.00, -1.16, 0.15, -1.25, -0.58, -0.09, 0.36, -0.08, 0.14, 0.32, -0.27, -1.06, -0.68, -0.85],
    [0.10, -0.09, -0.39, -0.47, 0.01, 0.19, 0.79, 0.15, 0.42, 0.13, 0.48, 0.26, 0.50, 0.15, 0.53, 0.10, -0.19, 0.10, 0.10, 0.04],
    [-0.36, -0.53, 0.28, 0.33, -1.08, 0.24, -0.24, -1.25, 0.13, -1.10, -0.50, 0.21, 0.42, -0.01, -0.07, 0.17, 0.07, -0.97, -0.95, -0.63],
    [-0.22, -0.43, 0.35, 0.29, -0.78, 0.02, -0.07, -0.58, 0.48, -0.50, -0.74, 0.32, 0.01, 0.26, 0.15, 0.48, 0.16, -0.73, -0.56, -1.02],
    [0.07, -0.52, -0.02, 0.01, 0.22, -0.04, 0.20, -0.09, 0.26, 0.21, 0.32, 0.14, 0.27, 0.37, 0.13, 0.15, 0.10, 0.40, -0.12, 0.32],
    [0.14, -0.14, 1.03, 0.69, 0.20, 0.60, 0.25, 0.36, 0.50, 0.42, 0.01, 0.27, 0.27, 1.02, 0.47, 0.54, 0.88, -0.02, -0.37, -0.12],
    [0.01, -0.43, 0.49, 0.04, 0.26, 0.46, 0.69, -0.08, 0.15, -0.01, 0.26, 0.37, 1.02, -0.12, 0.24, 0.29, 0.04, -0.11, 0.18, 0.11],
    [0.20, -0.24, -0.37, -0.52, -0.19, 0.50, 0.24, 0.14, 0.53, -0.07, 0.15, 0.13, 0.47, 0.24, 0.17, 0.27, 0.45, 0.01, -0.73, 0.01],
    [-0.09, 0.13, 0.19, 0.18, -0.22, 0.28, 0.21, 0.32, 0.10, 0.17, 0.48, 0.15, 0.54, 0.29, 0.27, -0.06, 0.08, 0.12, -0.22, -0.14],
    [-0.05, -0.22, -0.12, 0.37, 0.02, 0.28, 0.11, -0.27, 0.19, 0.07, 0.16, 0.10, 0.88, 0.04, 0.45, 0.08, -0.03, -0.01, 0.11, -0.32],
    [-0.42, -0.62, 0.69, 0.39, -1.15, 0.27, 0.16, -1.06, 0.10, -0.97, -0.73, 0.40, -0.02, -0.11, 0.01, 0.12, -0.01, -0.89, -0.56, -0.71],
    [0.05, 0.24, 0.04, 0.03, -0.60, 0.51, -0.85, -0.68, 0.10, -0.95, -0.56, -0.12, -0.37, 0.18, -0.73, -0.22, 0.11, -0.56, -0.05, -1.41],
    [-0.50, -0.79, 0.43, 0.17, -0.88, -0.35, -0.26, -0.85, 0.04, -0.63, -1.02, 0.32, -0.12, 0.11, 0.01, -0.14, -0.32, -0.71, -1.41, -0.76]
]
p_matrix = [
    [-1.65, -2.83, 1.16, 1.80, -3.73, -0.41, 1.90, -3.69, 0.49, -3.01, -2.08, 0.66, 1.54, 1.20, 0.98, -0.08, 0.46, -2.31, 0.32, -4.62],
    [-2.83, -39.58, -0.82, -0.53, -3.07, -2.96, -4.98, 0.34, -1.38, -2.15, 1.43, -4.18, -2.13, -2.91, -0.41, -2.33, -1.84, -0.16, 4.26, -4.46],
    [1.16, -0.82, 0.84, 1.97, -0.92, 0.88, -1.07, 0.68, -1.93, 0.23, 0.61, 0.32, 3.31, 2.67, -2.02, 0.91, -0.65, 0.94, -0.71, 0.90],
    [1.80, -0.53, 1.97, 1.45, 0.94, 1.31, 0.61, 1.30, -2.51, 1.14, 2.53, 0.20, 1.44, 0.10, -3.13, 0.81, 1.54, 0.12, -1.07, 1.29],
    [-3.73, -3.07, -0.92, 0.94, -11.25, 0.35, -3.57, -5.88, -0.82, -8.59, -5.34, 0.73, 0.32, 0.77, -0.40, -2.22, 0.11, -7.05, -7.09, -8.80],
    [-0.41, -2.96, 0.88, 1.31, 0.35, -0.20, 1.09, -0.65, -0.16, -0.55, -0.52, -0.32, 2.25, 1.11, 0.84, 0.71, 0.59, -0.38, 1.69, -1.90],
    [1.90, -4.98, -1.07, 0.61, -3.57, 1.09, 1.97, -0.71, 2.89, -0.86, -0.75, 1.84, 0.35, 2.64, 2.05, 0.82, -0.01, 0.27, -7.58, -3.20],
    [-3.69, 0.34, 0.68, 1.30, -5.88, -0.65, -0.71, -6.74, -0.01, -9.01, -3.62, -0.07, 0.12, -0.18, 0.19, -0.15, 0.63, -6.54, -3.78, -5.26],
    [0.49, -1.38, -1.93, -2.51, -0.82, -0.16, 2.89, -0.01, 1.24, 0.49, 1.61, 1.12, 0.51, 0.43, 2.34, 0.19, -1.11, 0.19, 0.02, -1.19],
    [-3.01, -2.15, 0.23, 1.14, -8.59, -0.55, -0.86, -9.01, 0.49, -6.37, -2.88, 0.97, 1.81, -0.58, -0.60, -0.41, 0.72, -5.43, -8.31, -4.90],
    [-2.08, 1.43, 0.61, 2.53, -5.34, -0.52, -0.75, -3.62, 1.61, -2.88, -6.49, 0.21, 0.75, 1.90, 2.09, 1.39, 0.63, -2.59, -6.88, -9.73],
    [0.66, -4.18, 0.32, 0.20, 0.73, -0.32, 1.84, -0.07, 1.12, 0.97, 0.21, 0.61, 1.15, 1.28, 1.08, 0.29, 0.46, 0.93, -0.74, 0.93],
    [1.54, -2.13, 3.31, 1.44, 0.32, 2.25, 0.35, 0.12, 0.51, 1.81, 0.75, 1.15, -0.42, 2.97, 1.06, 1.12, 1.65, 0.38, -2.06, -2.09],
    [1.20, -2.91, 2.67, 0.10, 0.77, 1.11, 2.64, -0.18, 0.43, -0.58, 1.90, 1.28, 2.97, -1.54, 0.91, 0.85, -0.07, -1.91, -0.76, 0.01],
    [0.98, -0.41, -2.02, -3.13, -0.40, 0.84, 2.05, 0.19, 2.34, -0.60, 2.09, 1.08, 1.06, 0.91, 0.21, 0.95, 0.98, 0.08, -5.89, 0.36],
    [-0.08, -2.33, 0.91, 0.81, -2.22, 0.71, 0.82, -0.15, 0.19, -0.41, 1.39, 0.29, 1.12, 0.85, 0.95, -0.48, -0.06, 0.13, -3.03, -0.82],
    [0.46, -1.84, -0.65, 1.54, 0.11, 0.59, -0.01, 0.63, -1.11, 0.72, 0.63, 0.46, 1.65, -0.07, 0.98, -0.06, -0.96, 1.14, -0.65, -0.37],
    [-2.31, -0.16, 0.94, 0.12, -7.05, -0.38, 0.27, -6.54, 0.19, -5.43, -2.59, 0.93, 0.38, -1.91, 0.08, 0.13, 1.14, -4.82, -2.13, -3.59],
    [0.32, 4.26, -0.71, -1.07, -7.09, 1.69, -7.58, -3.78, 0.02, -8.31, -6.88, -0.74, -2.06, -0.76, -5.89, -3.03, -0.65, -2.13, -1.73, -12.39],
    [-0.37, -4.62, -4.46, 0.90, 1.29, -8.80, -1.90, -3.20, -5.26, -1.19, -4.90, -9.73, 0.93, -2.09, 0.01, 0.36, -0.82, -3.59, -12.39, -2.68]
]
aa_list = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
p_matrix = np.array(p_matrix)
m_matrix = np.array(m_matrix)

def iupred(seq, sequence_separation=2, window_size=100, window_size_smooth=10):
    '''
    Calculate residue IUPRED energy considering neighbouring residues (windows_size) and
    smoothing by window_size_smooth
    :param seq: a string of aminoacids
    :param sequence_separation: neighbours min distance
    :param window_size: neighbours max distance
    :param window_size_smooth: sliding average window size
    :return: row prediction, smoothed prediction
    '''

    pred = []
    pred_smooth = []

    indices = [aa_list.index(aa) for aa in list(seq)]  # Transform sequence into indexes as in the P matrix
    for i, aa_index in enumerate(indices):

        # Get the slice i-100/i+100 excluding adjacent positions (+/-1)
        start_before = max(0, i - window_size)
        end_before = max(0, i - sequence_separation)
        start_after = min(len(indices) - 1, i + sequence_separation)
        end_after = min(len(indices) - 1, i + window_size)
        indices_local = indices[start_before: end_before] + indices[start_after: end_after]
        # print(i, aa_index, aa_list[aa_index], len(indices), len(indices_local), i, start_before, end_before, start_after, end_after)

        # Count amino acids in the window
        row = np.full((20,), 0)
        for index in indices_local:
            row[index] += 1
        # print(row)

        row = row / len(indices_local)  # calculate AA frequency
        # print(row)
        print(p_matrix[aa_index])
        row = row * p_matrix[aa_index, ]  # calculate energy
        # print(row)

        aa_energy = np.sum(row)
        # print(i, seq[i], aa_energy)

        pred.append(aa_energy)

    # Smooth the prediction (moving average)
    for i in range(len(pred)):
        frag = pred[max(0, i - window_size_smooth): min(i + window_size_smooth, len(pred))]
        pred_smooth.append(sum(frag) / len(frag))

    return pred, pred_smooth

def iupred2(seq, sequence_separation=2, window_size=100, window_size_smooth=10, matrix=[]):
    '''
    Calculate residue IUPRED energy considering neighbouring residues (windows_size) and
    smoothing by window_size_smooth
    :param seq: a string of aminoacids
    :param sequence_separation: neighbours min distance
    :param window_size: neighbours max distance
    :param window_size_smooth: sliding average window size
    :return: row prediction, smoothed prediction
    '''

    pred = []
    pred_smooth = []

    indices = [aa_list.index(aa) for aa in list(seq)]  # Transform sequence into indexes as in the P matrix

    for i, aa_index in enumerate(indices):

        # Get the slice i-100/i+100 excluding adjacent positions (+/-1)
        start_before = max(0, i - window_size)
        end_before = max(0, i - sequence_separation)
        start_after = min(len(indices) - 1, i + sequence_separation)
        end_after = min(len(indices) - 1, i + window_size)
        indices_local = indices[start_before: end_before] + indices[start_after: end_after]
        # print(i, aa_index, aa_list[aa_index], len(indices), len(indices_local), i, start_before, end_before, start_after, end_after)

        # Count amino acids in the window
        row = np.full((20,), 0)
        for index in indices:
            row[index] += matrix[index]

        #row = row / len(indices_local)  # calculate AA frequency
        # print(row)
        row = row * m_matrix[aa_index, ]  # calculate energy
        # print(row)

        aa_energy = np.sum(row)
        # print(i, seq[i], aa_energy)

        pred.append(aa_energy)

    # Smooth the prediction (moving average)
    for i in range(len(pred)):
        frag = pred[max(0, i - window_size_smooth): min(i + window_size_smooth, len(pred))]
        pred_smooth.append(sum(frag) / len(frag))

    return pred, pred_smooth

def question1() :
    
    pdbl = PDBList()
    pdbl.retrieve_pdb_file(pdb_id, pdir=path, file_format='pdb')
    structure = PDBParser(QUIET=True).get_structure(pdb_id, path + "pdb{}.ent".format(pdb_id))

    ns = NeighborSearch([atom for residue in structure[0]['A'] for atom in residue.get_atoms()])
    row = []
    for residue1, residue2 in ns.search_all(3.5, level="R"):  # level="R" returns pairs of residues in contact considering all atoms
        if residue1.id[0] == " " and residue2.id[0] == " ":  # Exclude hetero/water residues
            if abs(residue1.id[1] - residue2.id[1]) > 2:  # Sequence separation > 2
                row.append([residue1.id[1], residue2.id[1]])

    matrix = []
    for _ in structure[0]['A'] :
        array = []
        for _ in structure[0]['A'] :
            array.append(10)
        matrix.append(array)
    
    for r1,r2 in row :
        matrix[r1-1][r2-1] = 1
        matrix[r2-1][r1-1] = 1

    matrix = np.array(matrix, dtype = float)

    # Plot contact map
    contact_map = (matrix[:] < 8).astype(float)  
    # Calculate the contact map based on a distace threshold 8 Angstrom
    fig, ax = plt.subplots(figsize=(120, 120))
    im = ax.imshow(contact_map)

    # Set ticks
    ax.xaxis.set_major_locator(MultipleLocator(10))
    ax.xaxis.set_minor_locator(AutoMinorLocator(10))
    ax.yaxis.set_major_locator(MultipleLocator(10))
    ax.yaxis.set_minor_locator(AutoMinorLocator(10))

    plt.savefig(path + 'ca_contacts_{}.png'.format(pdb_id), bbox_inches='tight')
    return [residue for residue in structure[0]['A'] if residue.id[0] == " "], matrix

def question2(residues, matrix) :
    matrixNew = []
    for residue in matrix :
        integer = 0
        for residue1 in residue :
            if residue1 == 1 :
                integer += 1
        matrixNew.append(integer)
    matrixNew = np.array(matrixNew)

    seq = "".join([seq1(residue.get_resname()) for residue in residues])
    pred, pred_smooth = iupred2(seq, 2, 21, 10, matrixNew)

    fig, ax = plt.subplots(figsize=(12, 6))
    ax.set_title("{}_{}".format(pdb_id, chain_id))
    ax.axhline()
    ax.plot(np.arange(len(seq)), pred, ls='--')
    ax.plot(np.arange(len(seq)), pred_smooth, ls='-')

    plt.tight_layout()  # Remove figure padding
    plt.savefig(path + 'iupredExact_{}_{}.png'.format(pdb_id, chain_id), bbox_inches='tight')
    return pred

def question3(residues) :
    seq = "".join([seq1(residue.get_resname()) for residue in residues])
    pred, pred_smooth = iupred(seq, 2, 21, 10)
    
    fig, ax = plt.subplots(figsize=(12, 6))
    ax.set_title("{}_{}".format(pdb_id, chain_id))
    ax.axhline()
    ax.plot(np.arange(len(seq)), pred, ls='--')
    ax.plot(np.arange(len(seq)), pred_smooth, ls='-')

    plt.tight_layout()  # Remove figure padding
    plt.savefig(path + 'iupredEstimated_{}_{}.png'.format(pdb_id, chain_id), bbox_inches='tight')
    return pred

def __main__() :
    residues , matrix = question1()
    pred1 = question2(residues, matrix)
    pred2 = question3(residues)
    raw_1 = 0
    for i in pred1 :
        if i >= 0 :
            raw_1 += 1
    print(raw_1, raw_1/len(pred1))
    raw_2 = 0
    for i in pred2 :
        if i >= 0 :
            raw_2 += 1
    print(raw_2, raw_2/len(pred2))
__main__()

