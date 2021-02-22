# Rosetta Kinase CM.

Here are the STEPS and SCRIPTS used for the comparative modeling pipeline approach,

## 1. Conformers generation and ligand alignment (modeling_script.py)
* First, the maximum number of conformers will be generated for a given SMILES using OMEGA. The output will be a single SDF file with the maximum number of conformers.
* Second, the single SDF file will be used for the alignment of a query (Molecule of interest) and database (in-house active kinase ligand template library) molecules using ROCS. The output will be each template aligned query molecule.
* Third, combining the report files of each template aligned query molecule into a single report file.
* Fourth, based on the single report file, the top 100 conformers for the given SMILES are selected.
* Fifth, SDF to PARAMSwill be generated using 100 conformers for Rosetta Minimization step. 
```
python modeling_script.py -f ~/Desktop/2W1C_A_L0C -omega omega2 -rocs rocs -temp_lig /Users/kiruba/Desktop/rosetta_kinase_cm/template_ligand_library
```
## 2. Sequence alignment and protein modeling (new_protein_modeling.py)
* First, target-template sequence alignment will be performed using in-house active kinase sequence template library (EMBOSS Needleman-Wunsch algorithm)
* Second, selection of top hit templates will be applied using sequence and ligand similarity approach (defined as Template Score). Based on this approach, the top 10 templates for the given sequence will be selected.
* Third, 10 predicted models of target protein will be performed using PyRosetta.
* Forth, using the top first model we concatenate the 100 conformers from ligand alignment and this results into an unrefined protein-ligand complex of 100 comparative models.
## 3. Minimization of protein-ligand complex (minimization.py)
* First, the input files for Rosettaminimization process will be generated for parallel computing.
* Second, once minimization finished, the energy for each model will be calculated similarly to the first step.
## 4. Analysis (analysis_1.py)
* Here, I will generate the table for 100 minimized structures. It contains the name and energy attributes of those models. Out of 100, the top 10 models will be selected using Rosetta energy values.
## 5. Complex modeling of remaining protein models (top_comp_prtn_lig_modeling.py) (from step 2, third point)
* Here, the PARAMS files of top 10 ligands will be taken from step 1, fifth point.
* Concatenation of protein-ligand complex (this will again result into 100 complex models) 
## 6. Minimization (top_comp_prtn_lig_modeling_minimization.py)
* The minimization process is same as step 3.
## 7. Analysis (analysis_2.py)
* The analysis process is same as step 4.
* The top 1 model will be reported as the best prediction. 
For the ROCS, EMBOSS, and PyRosetta, I use parallel computing so I tried to keep those scripts separately.