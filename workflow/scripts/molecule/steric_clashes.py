#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Summarising the steric clashes in a pdb file.


"""
from Bio import PDB
import os
import sys 
import itertools
import logging
import polars as pl
from Bio.PDB import PDBParser 
import seaborn as sns
import numpy as np



# Configure the logging system
logging.basicConfig(level=logging.DEBUG,
                    format='%(asctime)s - %(levelname)s - %(message)s',
                    datefmt='%Y-%m-%d %H:%M:%S',
                    filename='app.log',
                    filemode='w')

console_handler = logging.StreamHandler(sys.stdout)
console_handler.setLevel(logging.DEBUG)
console_handler.setFormatter(
    logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s'))


def load_pdb(path: str):
    """Load a PDB structure

    Parameters
    ==========
    path:
        The path leading to the structure.
    """
    parser = PDBParser()
    structure = parser.get_structure(path,
                                     path)
    return structure


# Atomic radii for various atom types.
# You can comment out the ones you don't care about or add new ones
atom_radii = {
#    "H": 1.20,  # Who cares about hydrogen??
    "C": 1.70,
    "N": 1.55,
    "O": 1.52,
    "S": 1.80,
    "F": 1.47,
    "P": 1.80,
    "CL": 1.75,
    "MG": 1.73,
}

def count_clashes(structure, clash_cutoff=0.63):
    """
    source -- This code is adapted from a post made by Brennan Abanades Kenyon 
        on [this website](https://www.blopig.com/blog/2023/05/checking-your-pdb
        -file-for-clashing-atoms/)

    """
    # Set what we count as a clash for each pair of atoms
    clash_cutoffs = {i + "_" + j: (clash_cutoff * (atom_radii[i] +
                                                   atom_radii[j])) for i in
                     atom_radii for j in atom_radii}

    # Extract atoms for which we have a radii
    atoms = [x for x in structure.get_atoms() if x.element in atom_radii]
    coords = np.array([a.coord for a in atoms], dtype="d")

    # Build a KDTree (speedy!!!)
    kdt = PDB.kdtrees.KDTree(coords)

    # Initialize a list to hold clashes
    clashes = []

    # Iterate through all atoms
    for atom_1 in atoms:
        # Find atoms that could be clashing
        kdt_search = kdt.search(np.array(atom_1.coord, dtype="d"),
                                max(clash_cutoffs.values()))

        # Get index and distance of potential clashes
        potential_clash = [(a.index, a.radius) for a in kdt_search]

        for ix, atom_distance in potential_clash:
            atom_2 = atoms[ix]

            # Exclude clashes from atoms in the same residue
            if atom_1.parent.id == atom_2.parent.id:
                continue

            # Exclude clashes from peptide bonds
            elif (atom_2.name == "C" and atom_1.name == "N") or (atom_2.name ==
                                                                 "N" and
                                                                 atom_1.name ==
                                                                 "C"):
                continue

            # Exclude clashes from disulphide bridges
            elif (atom_2.name == "SG" and atom_1.name == "SG") and atom_distance > 1.88:
                continue

            if atom_distance < clash_cutoffs[atom_2.element + "_" + atom_1.element]:
                clashes.append((atom_1, atom_2))

    return len(clashes) // 2

def summarise_steric_clashes(structure):
    df_list = []
    for frame in structure.get_models():
        total_length =len([a for a in frame.get_atoms()])
        number = count_clashes(frame, clash_cutoff=0.63)
        df = (pl.DataFrame({"structure" : structure.id, "clashes" : number,
                         "clashes_percentage" : (number/total_length)*100}).
              with_columns(pl.lit(frame.id).alias("frameid"))
              )
        df_list.append(df)
    return pl.concat(df_list)

def main():
    # Get the root logger
    input_structure = snakemake.input[0]
    output_analysis = snakemake.output[0]
    logger = logging.getLogger()
    logger.addHandler(console_handler)
    path = input_structure
    pdb = load_pdb(path)
    df = summarise_steric_clashes(pdb)
    print(df)
    df.write_csv(output_analysis)

if __name__ == "__main__":

    main()

