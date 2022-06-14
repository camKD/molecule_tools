from rdkit import Chem
from rdkit.Chem.Lipinski import RotatableBondSmarts
from rdkit.Chem import rdchem # GetNeighbors
import numpy as np

def identify_all_rotatable_bonds(structure_path):
    mol = Chem.rdmolfiles.MolFromPDBFile(structure_path,removeHs=False)
    return identify_main_rotatable_bonds(mol) + identify_OH_rotatable_bonds(mol)

### either load mol directly or from path ###

def identify_main_rotatable_bonds(structure_path_or_rdkit_mol):
    dihedral_indexes = []
    mol = structure_path_or_rdkit_mol
    # if mol is passed at a pdb path, convert to mol. Otherwise use passed mol. 
    if isinstance(mol, str):
        mol = Chem.rdmolfiles.MolFromPDBFile(mol,removeHs=False)
    rotatableBonds = Chem.MolFromSmarts('[!$(*#*)&!D1]-&!@[!$(*#*)&!D1]')
    index_of_main_rotatable_bonds = mol.GetSubstructMatches(rotatableBonds)
    atoms = mol.GetAtoms()
    # Get neighbours
    for idx_pair in index_of_main_rotatable_bonds:
        # Get the two atoms
        atom_1 = atoms[idx_pair[0]]
        atom_2 = atoms[idx_pair[1]]
        # Get neighbours
        neighbours_1 = atom_1.GetNeighbors()
        neighbours_2 = atom_2.GetNeighbors()
        # Remove atom_1 from neighbours_2 and vice versa
        neighbours_1 = [atom for atom in neighbours_1 if atom.GetIdx() != atom_2.GetIdx()]
        neighbours_2 = [atom for atom in neighbours_2 if atom.GetIdx() != atom_1.GetIdx()]
        # Sort by heavy atoms, I think this ensures the least movement?
        neighbours_1 = sorted(neighbours_1, key=lambda atom: atom.GetMass(),reverse=True)
        neighbours_2 = sorted(neighbours_2, key=lambda atom: atom.GetMass(),reverse=True)
        # Take the heaviest atom
        atom_1_neighbour = neighbours_1[0]
        atom_2_neighbour = neighbours_2[0]
        # Add dihedral definition to list
        dihedral_indexes.append([atom_1_neighbour.GetIdx()]+list(idx_pair)+[atom_2_neighbour.GetIdx()])
    return dihedral_indexes
    

def identify_OH_rotatable_bonds(structure_path_or_rdkit_mol):
    dihedral_indexes = []
    mol = structure_path_or_rdkit_mol
    # if mol is passed at a pdb path, convert to mol. Otherwise use passed mol. 
    if isinstance(mol, str):
        mol = Chem.rdmolfiles.MolFromPDBFile(mol,removeHs=False)
    # add explicit hydrogens
    #mol = Chem.rdmolops.AddHs(mol)
    # Get index of Os in hydroxyls
    hydroxyls = Chem.MolFromSmarts('[OX2H]')
    index_of_hydroxyls = mol.GetSubstructMatches(hydroxyls)
    #print(index_of_hydroxyls)
    atoms = mol.GetAtoms()
    # iterate through Os
    for idx in index_of_hydroxyls:
        # get the oxygen atom
        O = atoms[idx[0]]
        # find neigbours ie the H and C opposite, could feasibly use this to get hydrogen angles in future
        neighbours = O.GetNeighbors()
        print(neighbours)
        #print(neighbours[0].GetIdx(),neighbours[1].GetIdx())
        neighbours = [atom for atom in neighbours]
        # Sort by light atoms so we start with H?
        neighbours = sorted(neighbours, key=lambda atom: atom.GetMass())
        print(neighbours)
        # Get neigbours of the first cardbon (to find more carbons) and remove the oxygen
        C_neighbours = neighbours[1].GetNeighbors()
        C_neighbours = [atom for atom in C_neighbours if atom.GetIdx() != O.GetIdx()]
        # Sort by mass
        C_neighbours = sorted(C_neighbours, key=lambda atom: atom.GetMass(),reverse=True)
        atom_1 = neighbours[0]
        atom_2 = O
        atom_3 = neighbours[1]
        atom_4 = C_neighbours[0]
        dihedral_indexes.append([atom_1.GetIdx(), atom_2.GetIdx(), atom_3.GetIdx(), atom_4.GetIdx()])
    return dihedral_indexes


### I'm an idiot. Just needed to add hydrogens to the original script. Ah well. identify_main_rotatable_bonds() will not identify all rotatable bonds

#print(identify_main_rotatable_bonds('epctcn.pdb'))
#print(identify_OH_rotatable_bonds('epctcn.pdb'))
#print(identify_all_rotatable_bonds('epctcn.pdb'))
