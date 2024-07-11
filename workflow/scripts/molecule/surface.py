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


def deterimine_sasa(path: str):
    """Get the SASA score per atom

    Parameters
    ==========
    path:
        The path leading to the structure.
    """
    protein = md.load(path)
    logger.info(f"loaded {protein} from {path}")
    sasa_data = md.shrake_rupley(protein, mode="residue")
    atoms, bonds = protein.top.to_dataframe()
    atoms = atoms[['chainID', 'resName', 'resSeq']].drop_duplicates()

    atoms.insert(1,"SASA", sasa_data[0], True)
    return pl.from_pandas(atoms)

def plot_sasa(df):
    return sns.lineplot(x='resSeq',
                 y='SASA',
                 data=df)


def write_sasa_plot(df, path):
    fig=plot_sasa(df).get_figure()
    fig.savefig(path)



def main():
    # Get the root logger
    input_structure = snakemake.input[0]
    output_plot = snakemake.output[0]
    output_analysis = snakemake.output[1]
    logger = logging.getLogger()
    logger.addHandler(console_handler)
    path = input_structure
    df = deterimine_sasa(path)
    write_sasa_plot(df, output_plot)
    df.write_csv(output_analysis)



if __name__ == "__main__":

    main()
