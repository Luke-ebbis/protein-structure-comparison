#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Code adapted from
import Bio.PDB
from Bio.SeqIO import PdbIO, FastaIO
from Bio import SeqIO
import os
import itertools
import argparse
from typing import List, Tuple, Set, Dict
import logging
import sys
logging.basicConfig(level=logging.DEBUG,
                    format='%(asctime)s - %(levelname)s - %(message)s',
                    datefmt='%Y-%m-%d %H:%M:%S',
                    filemode='w')
console_handler = logging.StreamHandler(sys.stdout)
console_handler.setLevel(logging.DEBUG)
console_handler.setFormatter(
logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s'))
logger = logging.getLogger()

def extract_sequences(pdb_filename,
                      name: str|None = None,
                      comment: str|None = None):
    if not name:
        name = os.path.basename(pdb_filename).replace('.pdb', "")
    if not comment:
        comment = "converted"
    sequences = []
    with open(pdb_filename, "r") as pdb_handle:
        seqs = SeqIO.parse(pdb_handle, "pdb-atom")
        for seq in seqs:
            sequences.append(seq)
    logger.info(f"Found {len(sequences)} chains in {pdb_filename}")
    fastas = [process_chain(s, name, comment) for s in sequences]
    return fastas

def process_chain(chain, name: str|None = None, comment: str|None = None):
    if chain.id.startswith("???"):
        chain_id = f"{name}_{chain.id.split(":")[1]}"
        chain.id = chain_id
    if chain.description.startswith("???"):
        chain_comment = comment
        chain.description = chain_comment
    return chain

def write_fasta(sequences, fasta_file):
    SeqIO.write(sequences, fasta_file, "fasta")
    logger.info(f"Wrote the sequences to {fasta_file}")

def main():
    parser = argparse.ArgumentParser(
                        prog='FASTA extractinator',
                        description=f'Extract the primary structure from a protein file',
                        epilog='Sibbe Bakker')
    parser.add_argument('structure',
                        help="PDB file.")
    parser.add_argument('output',
                        help="Name of the fastafile")
    args = parser.parse_args()
    pdb_filename = args.structure
    
    sequences = extract_sequences(args.structure)
    write_fasta(sequences, args.output)


if __name__ == "__main__":
    main()
