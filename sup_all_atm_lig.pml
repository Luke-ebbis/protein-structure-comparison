#!/usr/bin/env pymol
cmd.load("sup.pdb", "structure1")
cmd.load("results/input/6DV2.pdb", "structure2")
hide all
set all_states, off
show cartoon, structure1
show cartoon, structure2
show stick, not polymer
show sphere, not polymer
color blue, structure1
color red, structure2
set ribbon_width, 6
set stick_radius, 0.3
set sphere_scale, 0.25
set ray_shadow, 0
bg_color white
set transparency=0.2
zoom polymer and ((structure1 and c. A) or (structure2 and c. A))

