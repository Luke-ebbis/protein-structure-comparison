include: "convertion.smk"
include: "io.smk"

## contact_map:
##    Calculate the contact map of a structure.
##
rule contact_map:
  conda: "../envs/contacts.yml"
#"https://raw.githubusercontent.com/Luke-ebbis/distance-map/main/env.yml"
  input:
    data="results/input/{identifier}.pdb",
    script="../scripts/external/raw.githubusercontent.com/Luke-ebbis/distance-map/v0.1.0/distance-map.py",
  output: directory("results/contact-map/{identifier}")
  shell:
    """
    python3  {input.script} {input.data} {output}
    """

