"""Snakemake file for the conversion of structural biology formats
"""

from snakemake.utils import min_version
min_version("6.4.1")

## fasta:
##    Determine the primary structure of a protein.
##
rule fasta:
  conda: '../envs/comparison.yml'
  input: "results/input/{identifier}.pdb"
  output: "results/fasta/{identifier}.fasta"
  params:
    script=workflow.source_path("../scripts/pdb2fasta.py")
  shell:
    """
    python3 {params.script} {input} {output}
    """

## convert_cif2pdb:
##    Using Bio python, read in a CIF file, and write it out to a PDB file.
##  Adapted from Spencer Bliven (https://gist.github.com/sbliven/b7cc2c530)
##
rule convert_cif2pdb:
  """Convert a structure using biopython
  """
  conda: '../envs/comparison.yml'
  input: "results/input/{identifier}.cif"
  output: "results/input/{identifier}.pdb"
  params:
    script=workflow.source_path("../scripts/convert-cif-to-pdb.py")
  message: f"Converting {input} to PDB."
  shell: "python3 {params.script} {input} {output}"


