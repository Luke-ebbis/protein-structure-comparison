#!/usr/bin/env python
# -*- coding: utf-8 -*-
#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
import sys 
import numpy as np
import mdtraj as md
import itertools
import logging
import polars as pl
from Bio.PDB import PDBParser 
import seaborn as sns

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

def deterimine_beta_factors(structure):
    """Get the SASA score per atom

    Parameters
    ==========
    path:
        The path leading to the structure.
    """
    df = pl.DataFrame(schema={"frame":int,
                              "residue" : str,
                              "residue_number" : int,
                              "betafactor" : int})
    rows = []
    for model in structure.get_models():
        frame = model.full_id
        print(frame)
        for residue in model.get_residues():
            residue_name = residue.resname
            segment_id = residue.segid
            for atom in residue.get_atoms():
                atom_id = atom.full_id
                beta = atom.bfactor
                rows += [{ "frame": frame[1], "residue" : residue_name,
                          "residue_number" : str(residue.id[1]), "betafactor" :
                          float(beta)}]
    df = pl.DataFrame(rows)
    return df


def plot_beta_factor(df):
    # https://stackoverflow.com/questions/64052747/plot-with-seaborn-not-show-all-in-x-axis
    ax = sns.lineplot(x='residue_number', y='betafactor', data=df)
    ax.set_xticklabels([])  
    return ax


def write_beta_factor_plot(df, path):

    fig=plot_beta_factor(df).get_figure()
    fig.savefig(path)



def main():
    # Get the root logger
    input_structure = snakemake.input[0]
    output_plot = snakemake.output[0]
    output_analysis = snakemake.output[1]
    logger = logging.getLogger()
    logger.addHandler(console_handler)
    path = input_structure
    pdb = load_pdb(path)
    df = deterimine_beta_factors(pdb)
    print(df)
    write_beta_factor_plot(df, output_plot)
    df.write_csv(output_analysis)



if __name__ == "__main__":

    main()
