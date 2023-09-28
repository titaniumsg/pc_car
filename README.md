## Structure-based prediction of cross-reactive peptides
 
Here, we provide documentation and code for the *de novo* prediction of cross-reactive peptides using the 10LH:PHOX2B/HLA-A*24:02 complex as an example. 

For a detailed description of the implementation, please refer to our manuscript: Sun, Florio, & Gupta *et al*., Structural principles of peptide-centric Chimeric Antigen Receptor recognition guide therapeutic expansion. *bioRxiv* [doi.org/10.1101/2023.05.24.542108](https://doi.org/10.1101/2023.05.24.542108)

![Alt text](denovo.png?raw=true "Image")

### System requirements:
- MacOS/Linux (tested on MacOS 12.6.3)
- Rosetta (tested on 2020.08)
- Anaconda (tested on 4.11.0)
   - Python 3.8+ (tested on 3.8.15)
   - platformdirs
   - PyMOL
   - tqdm

### Installation:
1. Create a conda environment using the provided YAML file using the following command: `conda env create -f environment.yml` This step takes roughly 5 minutes on a normal desktop computer, depending on internet speed.
2. Obtain Rosetta via https://www.rosettacommons.org/software/license-and-download
3. (Optional) Download the latest version of HLA3DB via https://hla3db.research.chop.edu/ and directly replace the `MHC_pdbs` folder. Note that this repository contains a static version of HLA3DB downloaded on April 28th, 2023.

### Usage

1. Activate the conda environment: `conda activate pc_car`
2. Run the prediction: `python predict_cr_peptides.py`. The script takes ~40 minutes to run on a normal desktop computer.
