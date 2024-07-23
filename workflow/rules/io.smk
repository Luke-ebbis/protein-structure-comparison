"""Snakemake file for managing input/output to the pipeline

The main need for this file, is to dynamically grab an identifier from the
`data/` folder, or the rcsb-pdb.
"""

from snakemake.utils import min_version
min_version("6.4.1")

ruleorder: download_rcsb > convert_cif2pdb > move

## unzip
##    If you place zip files from alphafold3 in the data folder, this rule
##    extracts the structures.
##
checkpoint unzip:
  input: "data/{identifier}.zip"
  output:
    d=directory("data/{identifier}-unzipped")
  shell: "mkdir {output} -p; unzip {input} -d {output.d}"


## download_rcsb:
##    If a input file is missing, it is download with this rule.
##
checkpoint download_rcsb:
  output:
    "results/downloads/{identifier}.{ext}"
  shell:
    "wget -O {output} 'https://files.rcsb.org/download/{wildcards.identifier}.{wildcards.ext}' --no-check-certificate"


rule f:
  """
  This rule makes sure that the content of the AF3 folders gets detected.
  """
  input: rules.unzip.output.d
  output: "data/{identifier}-unzipped/{name}.cif"


def get_input(wildcards):
  """Determine the source of the input file

  branches:
    1) If the input file matches the name of a zip folder, then the
  folder gets unzipped, and the cif file gets extracted--assuming AF3 directory
  formats.
    2) File is present directly in the data folder. It gets moved directly to
  the input folder.
    3) The file needs to be downloaded.
  """
  # we find the file in the data folder
  final_target = wildcards.identifier
  final_format = wildcards.ext
  files = os.listdir("data")
  match = [f for f in files if wildcards.identifier.startswith(f.split(".")[0])]
  match = [f for f in files if f != ".gitkeep"]

  # If the model could be from AF3
  if "_model_" in wildcards.identifier:
    match = [f for f in match if f.startswith(wildcards.identifier.split("_model_")[0])]
  else:
    match = [f for f in match if f.startswith(wildcards.identifier)]
  if match:
    # file is local
    match = match.pop()
    if match.endswith("zip") and len(wildcards.identifier) != 4:
      unzip_ckp = checkpoints.unzip.get(identifier=match.replace(".zip", ""))
      print(unzip_ckp.output)
      final_target = f"{unzip_ckp.output[0]}/{wildcards.identifier}.{wildcards.ext}"
      return final_target
    else:
      final_target = f"data/{wildcards.identifier}"
    output = f"{final_target}.{final_format}"
  else:
    # We need to download the file
    download_rcsb_ckp = checkpoints.download_rcsb.get(**wildcards)
    output = download_rcsb_ckp.output
  return output


rule move:
  """
  This rule moves all data to the same folder, so the analysis
  rules are easier to write.
  """
  input: get_input
  output: "results/input/{identifier}.{ext}"
  shell: "cp {input} {output}"

rule get_external_script:
  output:
    "../scripts/external/{url}"
  shell:
    """
    wget -O {output}   "https://{wildcards.url}"
    """

rule get_external_env:
  output:
    "../envs/external/{url}"
  shell:
    """
    wget -O {output}   "https://{wildcards.url}"
    """

