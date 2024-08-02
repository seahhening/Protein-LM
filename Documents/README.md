# Analysis of Protein Sequence (for Generative Protein Design)

V0.2 by Xiang Yang, 31 July 2024

# Overview
AlphaFold3 uses a generative diffusion model to predict protein structure from input protein sequences. <br>

Protein designers can insert their own constraints to bias the diffusion model's gaussian denoising process, thereby producing novel proteins with desired properties. <br>

This project analyses a protein sequence on the basis of 5 parameters: N-Neighbour Count, Conservation Score, Hydrophobicity, Volume and Proton Donating/Accepting Ability.


# How to Use
For the selected protein, obtain the consensus sequence from the [AlphaFold Protein Database](https://alphafold.ebi.ac.uk/). This sequence will be used as a backbone for other experimentally determined PDB data to be mapped on.

Experimental information from PDB sequences are obtained from [Uniprot Database](https://www.uniprot.org/).

Save the pdb files in a folder named <i>protein_pdb</i> under your working directory. At the `main` function, update the variables <i>pdb_file_paths</i>, <i>ref_seq</i>, <i>protein_id</i> and <i>num_residues</i>.

Run the `main` function
<br>

## Example Use (Protein: Lactate Oxidase)

Uniprot PDB Files:
<img src="image.png" alt="Uniprot PDB Files" width=80%/>

PDB File Download from EMBL-EBI Page
<img src="image-2.png" alt="PDB Download_2du2" width=80%/>

AlphaFold Protein Database:
<img src="image-1.png" alt="AlphaFold Consensus Sequence" width=80%/>

#### Sample Inputs:
<img src="image-3.png" alt="Sample Inputs" width=80%/>

#### Sample Outputs
(1) Neighbour Count & Conservation Score Merged Plot
<img src="Documents/Compiled (Submission)/LOX_outputs/LOX_merged.png" alt="merged_LOX" width=80%/>

(2) Summary Plot
![LOX_summary_plot](Documents/Compiled (Submission)/LOX_outputs/LOX_summary.png)

(3) Summary in CSV
![alt text](image-4.png)

(4) Analysis Output in CSV


## Inputs: 
(1) Experimental PDB files of <i>target protein</i> obtained through Uniprot <br>
(2) Ref_Seq PDB file obtained through AlphaFoldDB <br>
(3) Number of amino acids <br>
(4) Weights for analysis of protein sequence

## Outputs:
(1) Compiled data of protein sequence with 5 parameters in csv format <br>
(2) Visual plot of the protein sequence and properties of amino acids <br>
(3) Weighted sums of normalised data for each amino acid


### Demo: (To be Updated)
Try out the [Notebook Demo](https://github.com/salesforce/LAVIS/blob/main/examples/blip2_instructed_generation.ipynb) of this function applied to DNMT1 (homo sapien)! [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/salesforce/LAVIS/blob/main/examples/blip2_instructed_generation.ipynb)

# Functions
### 1.1 download_conservation_scores
```
Downloads conservation score from Consurf Database using specific PDB file

Input:
- pdb_file_path (str): Path to PDB file
- protein_id (str): Protein ID
- chains (str): Chain ID

Output:
- Folder with conservation score in txt format, file named with PDB ID 
```
### 1.2 extract_consevation_score
```
Extracts conservation scores from downloaded conservation score files
```
### 1.3 extract_pdb_info
```
Reads a PDB file and extracts 3D spatial coordinates and amino acids.

Input:
- pdb_file_path (str): Path to the PDB file

Output:  
- A list of dictionaries, each containing information about each atom
```
### 2.1 find_nearest_neighbour
```
Find the n-neighbors closest to the 'CA' atoms in the data using brute force with Euclidean distance.

Input:
- data: list of dictionaries, each containing atom information
- angstrom (float): distance threshold for neighbour counting

Output:
- data: list of dictionaries, containing only 'CA' atoms with updated 'neighbour_count' field, and without 'model_id', 'atom_name', and 'atom_coords' keys.
```
### 2.2 find_hydrophobicity
```
Assign hydrophobicity values to an amino acid, using IMGT Scale

Input:
- aa_residue: Dict Object with 'residue_name'

Ouput:
- aa_residue: Dict Object with 'hydrophobicity' key
```
### 2.3 find_volume
```
Assign volume to an amino acid, using IMGT Scale

Input:
- aa_residue: Dict object with 'residue_name'

Ouput:
- aa_residue: Dict object with 'volume' key
```
### 2.4 find_proton_donor
```
Classifies amino acid based on proton accepting/donating ability

Input:
- aa_residue: Dict object with 'residue_name'

Ouput:
- aa_residue: Dict object with 'proton_donor' key
```
### 2.4 merge_residue_data
```
Merges conservation score and neighbour count, then assigns hydrophobicity, volume and proton donor/acceptor

Input:
- conservation_scores: List of dictionaries with conservation score
- neighbour_counts: List of dictionaries with neighbour count

Output:
- merged_data: List of dictionaries with residue id, residue name, consevation score, neighbour count, hydrophobicity, volume and proton donor/acceptor
```
### 2.5 pdb_to_compiled_vector
```
For each pdb sequence, outputs the combined data with all 5 parameters

Input:
- pdb_file_path (str)
- angstrom (int)
- protein_id (str)

Output:
- merged_residues: List of dictionaries with with residue id, residue name, consevation score, neighbour count, hydrophobicity, volume and proton donor/acceptor
```
### 2.6 merge_to_ref_seq
```
Updates cosesnsus sequence with information extracted from PDB file

Input:
- ref_seq (str): File path to consensus sequence
- completed (dict): Processed pdb file to be updates to ref_seq

Output:
- updated_ref_seq: List of dictionaries
```
### 3.1 plot_and_calculate_correlation
```
Plots neighbout count and conservation scores for all residues

Input: 
- data: List of dictionaries
- angstrom
- protein_id
- num_residues

Output:
- scatter plot: y axes are neighbour count and conservation score, x axis is residue id
- caculated correlation between neighbour counts and conservations scores
```
### 3.2 plot_residue_properties
```
Plots the different parameters of each amino acid residue

Input:
- data: List of dictionaries
- protein_id
- num_residues

Output:
- summary chart with 5 parameters
```
### 3.3 analyse_protein

### 3.4 analyse_protein2