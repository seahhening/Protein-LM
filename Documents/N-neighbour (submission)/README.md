# Amino N-neighbour Function (for Generative Protein Design)
<img src="https://github.com/Reishuen/Protein-LM/assets/102805134/168e1979-f202-4ae0-9305-d6ad24f7abed" width="50%"/>

V1.0 by Rei Shuen Ng, 5 July 2024

# Overview
AlphaFold3 uses a generative diffusion model to predict protein structure from input protein sequences. <br><br>
Protein designers can insert their own constraints to bias the diffusion model's gaussian denoising process, thereby producing novel proteins with desired properties. 

<div style="display: flex; justify-content: space-between;">
  <img src="https://github.com/Reishuen/Protein-LM/assets/102805134/842ea596-cedf-43d2-b664-846ef1c74bff" width="45%">
  <img src="https://github.com/Reishuen/Protein-LM/assets/102805134/f1d369aa-f7c3-4334-9201-f8d09400324d" width="45%">
</div>


# How to Use
The Amino N-neighbour Function's purpose is to calculate the distance between Carbon-Alpha atoms of amino acids within a protein, giving a N-neighbour count of each amino acid. Amino acids with the highest N-neighbour counts can be interpreted physically as having the most interactions with neighbouring amino acids, suggesting its criticality to the protein's function. <br>

Hence when designing novel proteins, we might want to conserve amino acids with N-neighbour counts.

This hypothesis is cross-validated against evolutionarily conserved regions of the protein, represented by conservation scores, which justify the biological importance of certain amino acids by showing its conservation across species.
<br><br>
<img src="https://github.com/Reishuen/Protein-LM/assets/102805134/684b1377-dbd2-4566-9a8d-987181fac204" width="45%"/>

# Example use
## Inputs: 
#### (1) Experimental PDB files of {target protein} obtained through Uniprot 
#### (2) RefSeq PDB file obtained through AlphaFoldDB

## Outputs:
#### (1) N-neighbour count per amino acid vector
#### (2) Distance matrix between CA atoms in amino acids
#### (3) Visualisation plots

<div style="display: flex; justify-content: space-between;">
  <img src="https://github.com/Reishuen/Protein-LM/assets/102805134/b76f6d5e-ccf9-4eb0-9faa-544569a71736" width="45%">
  <img src="https://github.com/Reishuen/Protein-LM/assets/102805134/229b0d3b-2cf6-453e-9540-a27b5e768be9" width="45%">
</div>


### Demo:
Try out the [Notebook Demo](https://github.com/salesforce/LAVIS/blob/main/examples/blip2_instructed_generation.ipynb) of this function applied to DNMT1 (homo sapien)! [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/salesforce/LAVIS/blob/main/examples/blip2_instructed_generation.ipynb)

# Functions
### 1.1) pdb_to_n_neighbour_vector
```
pdb_to_n_neighbour_vector(pdb_file_path, angstrom = 5, save_distance_matrix = True)
```
### 1.2) merge_to_ref_seq
```
# Call this merging function after processing each Experimental PDB seq
merge_to_ref_seq(ref_seq, other_seq)
```

### 1.3) visualise_seq (TBC)
```
visualise_seq(ref_seq, count >= 5)
```

### 1.4) plot_conservation_vs_neighbours (for cross-validation)
```
plot_conservation_vs_neighbours(updated_ref_seq)
```

