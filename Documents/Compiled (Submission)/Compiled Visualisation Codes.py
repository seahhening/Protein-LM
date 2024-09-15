import pymol
from pathlib import Path
import json
import requests
import os
import pandas as pd
import math
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import numpy as np
import re

#input area: edit the information about the protein of interest 
protein_id = 'LOX' 
uniprot_id = 'Q44467'
ref_seq_directory = f'Documents/Compiled (Submission)/{protein_id}_pdb/AF-{uniprot_id}-F1-model_v4.pdb'
ref_DSSP_filepath = f'Documents/Compiled (Submission)/{protein_id}_dssp/AF-{uniprot_id}-F1-model_v4.dssp'
summary_csv_filepath = f'Documents/Compiled (Submission)/{protein_id}_outputs/{protein_id}_summary.csv'

# input the URL to download the json file
url =f'https://rest.uniprot.org/uniprotkb/{uniprot_id}.json?fields=ft_binding%2Cft_site%2Cft_act_site'

#Edit the switches to decide the type of output you want!
annotate_sites_pse = True
annotate_neighbours_pse = True
annotate_hydrophobicity_pse = True
annotate_dssp_pse = True

#Output 1: PyMOL session file to annotate active sites, binding sites, residues with high neighbour count and their respective neighbouring residues
ref_json_directory = f'Documents/Compiled (Submission)/{protein_id}_pdb/{uniprot_id}.json'

def sanitize_name(ligand_name):
    sanitized_ligand_name = re.sub(r'[^\w]', '_', ligand_name)
    return sanitized_ligand_name

#create a dictionary of lists to generate the input required for PyMOL to recognise the sites
with open(ref_json_directory, 'r') as file:
        data =json.load(file)
        active_site_residues = []
        binding_site_residues = {}
        unclassified_residues = []
        for feature in data.get('features', []):
            feature_type = feature.get('type')
            start = feature['location']['start']['value']
            end = feature['location']['end']['value']
            residue_range = list(range(start, end + 1))

            if feature_type == "Active site":
                active_site_residues.extend(residue_range)
            elif feature_type == "Binding site":
                ligand_name = feature.get('ligand', {}).get('name', 'Unknown_Ligand')
                
                # Sanitize ligand name to make it PyMOL-compatible
                sanitized_ligand_name = sanitize_name(ligand_name)

                if sanitized_ligand_name not in binding_site_residues:
                    binding_site_residues[sanitized_ligand_name] = []
                binding_site_residues[sanitized_ligand_name].extend(residue_range)
            elif feature_type == "Site":
                unclassified_residues.extend(residue_range)
            else:
                print("No relevant sites have been identified.")

select = {'ACTIVE SITE RESIDUES': active_site_residues,
'UNCLASSIFIED RESIDUES' : unclassified_residues}

selection_dict ={}

#stores each binding site as separate key-value pair entries in the dictionary 
for ligand_name, residues in binding_site_residues.items():
    selection_dict[f'BINDING SITE RESIDUES {ligand_name}'] = residues

selection_dict.update(select)


#generate strings that can be PyMOL to annotate the sites 
for site_type, residues in selection_dict.items():
    selection_dict[site_type] = "+".join(map(str, residues))
print(selection_dict)

#function to annotate the sites
def annotate_sites(selection_dict, annotate_sites_pse=True):
    if annotate_sites_pse:
        for site_type, selection_str in selection_dict.items():
            try:
                # Create a unique selection name for each site type
                site_name = site_type.replace(' ', '_').lower()
                
                # Ensure the previous color is cleared
                pymol.cmd.delete(site_name)

                # Select residues and show them
                pymol.cmd.select(site_name, f'resi {selection_str}')
                pymol.cmd.show('sticks', site_name)

                # Color and label residues based on site type
                if 'binding_site_residues' in site_name:
                    pymol.cmd.color('pink', site_name)  # Same color for all binding sites
                elif 'active_site_residues' in site_name:
                    pymol.cmd.color('hotpink', site_name)
                elif 'unclassified_residues' in site_name:
                    pymol.cmd.color('dirtyviolet', site_name)

                pymol.cmd.label(f"byres ({site_name})", "resn + resi")
                pymol.cmd.set("label_color", "black", site_name)
                pymol.cmd.set("label_size", -0.5)

            except Exception as e:
                print(f"An error occurred with selection of {site_type}: {e}")
    else:
        print('Annotation of active sites and binding sites is switched off.')

#open csv file to obtain the residues with neighbour count>5
my_residues = []
with open(summary_csv_filepath, 'r') as csv_file:
    df = pd.read_csv(summary_csv_filepath)
    filtered_df = df[df['neighbour_count']>=5] #change the threhsold if required
    my_residues.extend(filtered_df['residue_id'])


#generate selections for alpha carbons of residues that are within 5Å of the selected residuesw
resi_selections = []
neighbors_CA_selections = []
def generate_neighbours(my_residues):
    for resi_num in my_residues:
        current_residue_CA_info = []
        neighbors_CA_info = []

        # Select the current residue's alpha carbon
        pymol.cmd.select(f"current_residue_{resi_num}_CA", f"name CA and resi {resi_num}")
        pymol.cmd.iterate(f"current_residue_{resi_num}_CA", "current_residue_CA_info.append((resi))", space = {'current_residue_CA_info':current_residue_CA_info})

        # Select only the Cα atoms of neighbors within 5Å of the current residue's alpha carbon
        pymol.cmd.select(f"neighbors_CA_{resi_num}", f"(name CA within 5.0 of current_residue_{resi_num}_CA) and not (resi {resi_num})")
        pymol.cmd.iterate(f"neighbors_CA_{resi_num}", "neighbors_CA_info.append((resi))", space= {'neighbors_CA_info':neighbors_CA_info})
        
        # Append to selection lists
        resi_selections.append(current_residue_CA_info)
        neighbors_CA_selections.append(neighbors_CA_info)
    return resi_selections, neighbors_CA_selections


def annotate_neighbours(neighbour_dict): 
    if annotate_neighbours_pse: 
        for resi_selection, neighbors_CA_string in neighbour_dict.items():
            try:
                # Select the current residue's alpha carbon and show it
                pymol.cmd.select(f"current_residue_{resi_selection}_CA", f"resi {resi_selection} and name CA")
                pymol.cmd.select(f"current_residue_{resi_selection}", f"resi {resi_selection}")
                pymol.cmd.show("sticks", f"current_residue_{resi_selection}")
                pymol.cmd.color("yellow", f"current_residue_{resi_selection}")
                pymol.cmd.color("green", f"current_residue_{resi_selection}_CA")
                
                pymol.cmd.label(f"current_residue_{resi_selection}_CA", 'resn')
                pymol.cmd.set('label_size', 12)  # Adjust label size as needed
                pymol.cmd.set('label_color', 'white')  # Adjust label color for visibility
                pymol.cmd.set('label_font_id', 7)
                
                # Process neighbors and select them
                neighbors_CA_list = neighbors_CA_string.split('+')
                for neighbor_resi in neighbors_CA_list:
                    pymol.cmd.select(f"neighbor_CA_{neighbor_resi}", f"resi {neighbor_resi} and name CA")
                    pymol.cmd.show("sticks", f"neighbor_CA_{neighbor_resi}")
                    pymol.cmd.color("cyan", f"neighbor_CA_{neighbor_resi}")
                    
                    pymol.cmd.distance(f"dist_{resi_selection}", f"current_residue_{resi_selection}_CA", f"neighbor_CA_{neighbor_resi}")
                
            except Exception as e:
                print(f"An error occurred with residue '{resi_selection}' and neighbors '{neighbors_CA_string}': {e}")
                return
    else:
        print('Annotation of neighbours is switched off. ')


def site_visualisation_session(ref_seq_directory, ref_json_directory, my_residues): 
    pymol.finish_launching(['pymol', '-cq']) #change this HEADLESS mode later
    pymol.cmd.load(ref_seq_directory)  # Load the PDB file

    # Set up the view and rendering options
    pymol.cmd.set("antialias", 1)
    pymol.cmd.set("orthoscopic", 1)
    pymol.cmd.set("gamma", 1.15)
    pymol.cmd.set("cartoon_fancy_helices", 1)
    pymol.cmd.set("cartoon_fancy_sheets", 1)
    pymol.cmd.color("gray80")
    pymol.cmd.set("ray_shadows", 0)
    pymol.cmd.set("ray_trace_fog", 1)

    #function to determine the nearest neighbours
    resi_selections = []
    neighbors_CA_selections = []
    resi_selections, neighbors_CA_selections = generate_neighbours(my_residues)

    neighbour_dict = {}
    for resi_selection, neighbors_CA_selection in zip(resi_selections, neighbors_CA_selections):
        key = ''.join(resi_selection)
        value = '+'.join(neighbors_CA_selection)
        neighbour_dict[key] = value #ensure the keys and values pair up
    
    #apply the function to annotate the neighbours
    annotate_neighbours(neighbour_dict=neighbour_dict)

    #apply the function to annotate sites
    annotate_sites(selection_dict = selection_dict)

    output_dir = Path(f'Documents/Compiled (Submission)/{protein_id}_outputs')
    output_dir.mkdir(parents=True, exist_ok=True)
    session_path = output_dir/(f'{protein_id}_site and neighbour visualisation.pse')
    pymol.cmd.save(str(session_path))
    pymol.cmd.do('reinitialize')


#Output 2: Visualisation of hydrophobic patches present on the protein
# Read and clean data from CSV
df2 = pd.read_csv(summary_csv_filepath)
resi_ids = df2['residue_id']
hydrophobicity_values = df2['hydrophobicity']

# Remove NaN values as no colour can be attributed
hydrophobicity_values_clean = [x for x in hydrophobicity_values if not math.isnan(x)]
hydrophobicity_dict = dict(zip(resi_ids, hydrophobicity_values_clean))
print(hydrophobicity_dict)


# Function to attribute hydrophobicity values to a color
def hydrophobicity_to_colour(hydrophobicity_value):
    hydrophobic_score = 4.5  # Max hydrophobic score based on IMGT scale
    hydrophilic_score = -4.5 # based on IMGT scale 

    cmap = plt.cm.rainbow

    # Normalize the value between 0 and 1
    normalized_value = (hydrophobicity_value - hydrophilic_score) / (hydrophobic_score - hydrophilic_score)
    normalized_value = np.clip(normalized_value, 0, 1)  # Ensure value is within [0, 1]

    # Get the color from the rainbow colormap
    color = cmap(normalized_value) #calls the normal normalised value to get a colour from the colourmap 

    # Convert from RGBA to RGB
    return [color[0], color[1], color[2]]

# Function to annotate protein structure with hydrophobicity coloring
def annotate_hydrophobicity(ref_seq_directory, hydrophobicity_dict):
    if annotate_hydrophobicity_pse:

        pymol.finish_launching(['pymol', '-cq'])  # Launch PyMOL in headless mode
        pymol.cmd.load(ref_seq_directory)  # Load the PDB file

        # Rendering options for better visualization
        pymol.cmd.set("antialias", 1)
        pymol.cmd.set("orthoscopic", 1)
        pymol.cmd.set("gamma", 1.15)
        pymol.cmd.set("cartoon_fancy_helices", 1)
        pymol.cmd.set("cartoon_fancy_sheets", 1)
        pymol.cmd.set("ray_shadows", 0)
        pymol.cmd.set("ray_trace_fog", 1)

        # Loop through each residue and apply coloring based on hydrophobicity
        for resi_id, hydrophobicity_value in hydrophobicity_dict.items():
            colour = hydrophobicity_to_colour(hydrophobicity_value)  # Get color based on hydrophobicity
            pymol.cmd.set_color(f"color_{resi_id}", colour)  # Define color in PyMOL
            pymol.cmd.select(f"resi_{resi_id}", f"resi {resi_id}")  # Select the residue by ID
            pymol.cmd.color(f"color_{resi_id}", f"resi_{resi_id}")  # Apply color to the residue
        
        # Save the PyMOL session to the output directory
        output_dir = Path(f'Documents/Compiled (Submission)/{protein_id}_outputs')
        output_dir.mkdir(parents=True, exist_ok=True)
        session_path = output_dir / f'{protein_id}_hydrophobicity_visualization.pse'
        pymol.cmd.save(str(session_path))
        pymol.cmd.do('reinitialize')
    else:
        print('Annotation of hydrophobicity is switched off. ')

#Output3 Visualisation of Proteins based on DSSP values
residue_number =[]
relative_solvent_accessibility = []

with open(ref_DSSP_filepath, 'r') as file:
    for line in file:
        # Skip lines that do not contain residue data (e.g., header or footer lines)
        if line.startswith(' '): 
            continue
        resi_number = line[9:16].strip()
        residue_number.append(resi_number)

        relative_ASA = line[42:49].strip()
        relative_solvent_accessibility.append(relative_ASA)


del(residue_number[0])
del(relative_solvent_accessibility[0])
relative_solvent_accessibility = [float(rel_SA) for rel_SA in relative_solvent_accessibility]
#the closer the value of REL ASA is to 1, the greater the solvent accessible solvent area
dssp_dict = dict(zip(residue_number, relative_solvent_accessibility ))


#function to map the DSSP values 
def dssp_to_colour(relative_solvent_accessibility):
    cmap = plt.cm.rainbow
    colour = cmap(relative_solvent_accessibility) #currently stored in RGBA format 
    return [colour[0], colour[1], colour[2]] #extracts RGB data from RGBA data

#function to open and annotate the PyMOL session file 
def annotate_dssp(ref_seq_directory, dssp_dict):
    if annotate_dssp_pse:

        pymol.finish_launching(['pymol', '-cq'])  # Launch PyMOL in headless mode
        pymol.cmd.load(ref_seq_directory)  # Load the PDB file

        # Rendering options for better visualization
        pymol.cmd.set("antialias", 1)
        pymol.cmd.set("orthoscopic", 1)
        pymol.cmd.set("gamma", 1.15)
        pymol.cmd.set("cartoon_fancy_helices", 1)
        pymol.cmd.set("cartoon_fancy_sheets", 1)
        pymol.cmd.set("ray_shadows", 0)
        pymol.cmd.set("ray_trace_fog", 1)

        for residue_number, relative_solvent_accessibility in dssp_dict.items():
            colour = dssp_to_colour(relative_solvent_accessibility = relative_solvent_accessibility)
            pymol.cmd.set_color(f"color_{residue_number}", colour)
            pymol.cmd.select(f"resi_{residue_number}", f"resi {residue_number}")  # Select the residue by ID
            pymol.cmd.color(f"color_{residue_number}", f"resi_{residue_number}")

        output_dir = Path(f'Documents/Compiled (Submission)/{protein_id}_outputs')
        output_dir.mkdir(parents=True, exist_ok=True)
        session_path = output_dir / f'{protein_id}_dssp_visualization.pse'
        pymol.cmd.save(str(session_path))
        pymol.cmd.do('reinitialize')
    else:
        print('Annotation of DSSP is switched off')

#main control function 
def generate_all_output_files(ref_seq_directory, ref_json_directory, my_residues,hydrophobicity_dict, dssp_dict):
    site_visualisation_session(ref_seq_directory=ref_seq_directory, ref_json_directory=ref_json_directory, my_residues=my_residues)
    annotate_hydrophobicity(ref_seq_directory=ref_seq_directory, hydrophobicity_dict=hydrophobicity_dict)
    annotate_dssp(ref_seq_directory = ref_seq_directory, dssp_dict = dssp_dict)

main_output = generate_all_output_files(ref_seq_directory=ref_seq_directory, ref_json_directory=ref_json_directory, my_residues=my_residues, hydrophobicity_dict=hydrophobicity_dict, dssp_dict = dssp_dict)