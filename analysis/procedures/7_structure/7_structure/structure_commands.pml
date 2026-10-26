# PyMOL commands, typed one by one in the PyMOL command line (PyMOL> prompt).
# Files 2VGB.pdb and 2VGI.pdb were downloaded from the RCSB PDB (PDB format; copies in input/) into the PyMOL working folder.

# Step 1. Load the wild-type (2VGB) and the R486W mutant (2VGI) structures
load 2VGB.pdb
load 2VGI.pdb
remove solvent
bg_color white
set ray_opaque_background, 1

# Step 2. Keep subunit A of each structure (both files contain four subunits, A-D)
remove not chain A
orient 2VGB
turn x, -90

# Step 3. Superpose the mutant on the wild type (all 517 C-alpha atoms, no atoms excluded)
align 2VGI and name CA, 2VGB and name CA, cycles=0

# Step 4. Colour the domains (residue limits from Valentini et al. 2002)
color grey70, all
color yellow, resi 57-84
color skyblue, resi 85-159+263-431
color orange, resi 160-262
color palegreen, resi 432-574
set cartoon_transparency, 0.6, 2VGI

# Step 5. Show residue 486, the substrate analogue phosphoglycolate (PGA) and FBP
show spheres, 2VGB and resi 486 and not name N+C+O
color red, 2VGB and resi 486 and elem C
show sticks, 2VGB and resn PGA+FBP
color black, 2VGB and resn PGA and elem C
color magenta, 2VGB and resn MN
show spheres, 2VGB and resn MN

# Step 6. Distance from Arg486 to the active site (phosphoglycolate)
distance d_active, 2VGB and resi 486 and name CZ, 2VGB and resn PGA and name C1
get_distance 2VGB and resi 486 and name CZ, 2VGB and resn PGA and name C1

# Step 7. Close-up of residue 486 and the A-domain loop (Asn361, Leu362)
hide spheres
show sticks, resi 486 and not name N+C+O
show sticks, resi 361+362
color red, 2VGB and resi 486 and elem C
color purple, 2VGI and resi 486 and elem C
disable d_active
set cartoon_transparency, 0.75
orient resi 486+361+362
zoom resi 486+361+362, 4.5

# Step 8. Hydrogen bonds of residue 486 with the 361-362 backbone
distance hb_wt, 2VGB and resi 486 and name NH2, 2VGB and resi 362 and name O
get_distance 2VGB and resi 486 and name NH2, 2VGB and resi 362 and name O
distance hb_mut1, 2VGI and resi 486 and name NE1, 2VGI and resi 361 and name O
get_distance 2VGI and resi 486 and name NE1, 2VGI and resi 361 and name O
distance hb_mut2, 2VGI and resi 486 and name NE1, 2VGI and resi 362 and name O
get_distance 2VGI and resi 486 and name NE1, 2VGI and resi 362 and name O

# Step 9. Figure 1: one subunit with its domains, residue 486 and the active site
set ray_opaque_background, 1
enable d_active
hide everything, 2VGI
hide sticks, resi 361+362
show spheres, 2VGB and resi 486 and not name N+C+O
color red, 2VGB and resi 486
show spheres, 2VGB and resn PGA+MN
color black, 2VGB and resn PGA
color magenta, 2VGB and resn MN
show sticks, 2VGB and resn FBP
color yellow, 2VGB and resn FBP
set cartoon_transparency, 0, 2VGB
orient 2VGB
turn x, -90
turn y, 180
zoom 2VGB, -2
set float_labels, on
set label_size, 26
set label_font_id, 7
set label_color, black
set dash_color, black, d_active
set dash_width, 4
set label_distance_digits, 0
set label_size, 34, d_active
set label_position, (0, 0, 20), d_active
set label_color, black, d_active
pseudoatom lbl_A, selection=2VGB and resi 85-159+263-431 and name CA
pseudoatom lbl_B, selection=2VGB and resi 160-262 and name CA
pseudoatom lbl_C, selection=2VGB and resi 432-574 and name CA
pseudoatom lbl_N, selection=2VGB and resi 57-84 and name CA
pseudoatom lbl_486, selection=2VGB and resi 486 and name CZ
pseudoatom lbl_PGA, selection=2VGB and resn PGA
pseudoatom lbl_FBP, selection=2VGB and resn FBP
hide everything, lbl_*
label lbl_A, "A domain"
label lbl_B, "B domain"
label lbl_C, "C domain"
label lbl_N, "N domain"
label lbl_486, "Arg486"
label lbl_PGA, "Active site"
label lbl_FBP, "FBP site"
set label_color, skyblue, lbl_A
set label_color, orange, lbl_B
set label_color, forest, lbl_C
set label_color, grey40, lbl_N
set label_color, red, lbl_486
set label_color, black, lbl_PGA
set label_color, olive, lbl_FBP
set label_position, (16, 54, 20), lbl_A
set label_position, (2, 24, 20), lbl_B
set label_position, (-34, 14, 20), lbl_C
set label_position, (-12, 30, 20), lbl_N
set label_position, (-2, -12, 20), lbl_486
set label_position, (10, -34, 20), lbl_PGA
set label_position, (-14, -18, 20), lbl_FBP
ray 2400, 1700
png Figure1_subunit_domains.png, dpi=300

# Step 10. Figure 2: residue 486 in the wild type and in the mutant, superposed
hide labels, lbl_A or lbl_B or lbl_C or lbl_N or lbl_486 or lbl_PGA or lbl_FBP
disable d_active
hide spheres
hide sticks
show cartoon, 2VGI
set cartoon_transparency, 0.8, 2VGB
set cartoon_transparency, 0.8, 2VGI
show sticks, resi 486 and not name N+C+O
show sticks, resi 361+362
set_color wt_red, [0.85, 0.0, 0.0]
set_color mut_purple, [0.45, 0.0, 0.6]
color wt_red, 2VGB and resi 486 and elem C
color mut_purple, 2VGI and resi 486 and elem C
color skyblue, resi 361+362 and elem C
color blue, (resi 361+362+486) and elem N
color red, (resi 361+362) and elem O
set label_distance_digits, 1
set dash_width, 5
set dash_color, wt_red, hb_wt
set dash_color, mut_purple, hb_mut1
set dash_color, mut_purple, hb_mut2
set label_color, wt_red, hb_wt
set label_color, mut_purple, hb_mut1
set label_color, mut_purple, hb_mut2
pseudoatom lbl_R, selection=2VGB and resi 486 and name CB
pseudoatom lbl_W, selection=2VGI and resi 486 and name CZ2
pseudoatom lbl_361, selection=2VGB and resi 361 and name CB
pseudoatom lbl_362, selection=2VGB and resi 362 and name CD1
hide everything, lbl_R or lbl_W or lbl_361 or lbl_362
label lbl_R, "Arg486 (wild type)"
label lbl_W, "Trp486 (R486W)"
label lbl_361, "Asn361"
label lbl_362, "Leu362"
set label_color, wt_red, lbl_R
set label_color, mut_purple, lbl_W
set label_color, black, lbl_361
set label_color, black, lbl_362
set label_position, (2, 7, 4), lbl_R
set label_position, (-7, -6, 4), lbl_W
set label_position, (5, 4, 4), lbl_361
set label_position, (6, -4, 4), lbl_362
orient resi 486+361+362
zoom resi 486+361+362, 6
set label_size, 22, hb_wt
set label_size, 22, hb_mut1
set label_size, 22, hb_mut2
ray 2400, 1700
png Figure2_R486W_closeup.png, dpi=300
