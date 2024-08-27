import pymol
from pymol import cmd

# Launch PyMOL
pymol.finish_launching()

# Initial display settings
cmd.set("antialias", 1)
cmd.set("orthoscopic", 1)
cmd.set("gamma", 1.15)
cmd.set("cartoon_fancy_helices", 1)
cmd.set("cartoon_fancy_sheets", 1)
cmd.set_color("wgray", [0.800, 0.800, 0.800])
cmd.set("ray_shadows", 0)
cmd.set("ray_trace_fog", 1)

# Color the whole structure by element in gray
cmd.color("wgray", "all")

# List of residues to analyze
my_residues = [28, 29, 30, 136, 142, 259, 327, 329]

# Initialize empty lists to hold selections for residues and neighbors
resi_selections = []
neighbors_CA_selections = []

# Loop through each residue
for resi_num in my_residues:
    # Select the current residue's alpha carbon
    cmd.select(f"current_residue_{resi_num}_CA", f"resi {resi_num} and name CA")
    
    # Select the whole residue
    cmd.select(f"current_residue_{resi_num}", f"resi {resi_num}")
    
    # Select only the Cα atoms of neighbors within 5Å of the current residue's alpha carbon
    cmd.select(f"neighbors_CA_{resi_num}", f"(name CA within 5.0 of current_residue_{resi_num}_CA) and not (resi {resi_num})")

    # Measure distances between the current residue's CA and its neighboring CAs
    cmd.distance(f"dist_{resi_num}", f"current_residue_{resi_num}_CA", f"neighbors_CA_{resi_num}")

    # Append to selection lists
    resi_selections.append(f"current_residue_{resi_num}")
    neighbors_CA_selections.append(f"neighbors_CA_{resi_num}")

# Combine selections into a single string
resi_selection = " or ".join(resi_selections)
neighbors_CA_selection = " or ".join(neighbors_CA_selections)

# Show the whole residues and their neighbors as sticks
cmd.show("sticks", resi_selection)
cmd.show("sticks", neighbors_CA_selection)

# Color neighbors' Cα cyan
cmd.color("cyan", f"({neighbors_CA_selection}) and name CA")

# Color Cα of residues green
cmd.color("green", f"({resi_selection}) and name CA")

# Color the whole residue yellow
cmd.color("yellow", f"({resi_selection}) and elem C")

# Debugging: Verify the selections
print("Residue Selection: ", resi_selection)
print("Neighbors CA Selection: ", neighbors_CA_selection)

# Label only the CA atoms of residues with their 3-digit code
cmd.label(f"({resi_selection}) and name CA", 'resn')

# Adjust label settings for better visibility
cmd.set('label_size', 12)  # Adjust label size as needed
cmd.set('label_color', 'white')  # Adjust label color for visibility
cmd.set('label_font_id', 7)  # Adjust font type if needed

# Update the view
cmd.zoom()
