# Mass Spectrometry Data Cleaning for Improved Identification of Post-Translational Modifications
## Data
### LC-MS data
For the development of the calibration module, .raw files from the following PRIDE projects was used and converted to mzML format:
- PXD001468 (https://www.ebi.ac.uk/pride/archive/projects/PXD001468)
- PXD032235 (https://www.ebi.ac.uk/pride/archive/projects/PXD032235)

### Search parameters
The following Sage search parameter JSON files were used and can be found in the "data" folder:
- Initial closed searches: SearchClosed.json
- Open search window selection: SearchClosed5ppm.json to SearchClosed30ppm.json
- Open searches: SearchOpen.json

### Sequence database
One FASTA file was consistently used through all the searches and can also be found in the "data" folder:
- human_ref_prot_and_contaminants.fasta

The custom Python script to merge the human and contaminant sequences can be found in the "scripts" folder:

## Mass calibration pipeline
All Python an Jupyter Notebook files written for the construction of the mass calibration module can also be found in the "scripts" folder.

The backbonde of the calibration pipeline itself was based on the 'Mass Spectrometry Workflow' Nextflow tool from GitHub user jvdnheme (https://github.ugent.be/jvdnheme/automatic_annotation). Extra modules including functionality to convert .xz into .raw files and to calibrate mzML files (LR and RF) were incorporated. The Nextflow script can be found in the "CalibScript.nf" file.

## Pipeline Requirements and Installation
...
## Pipeline Schematic Representation
![pipeline](https://github.com/NathanMRX2404/Thesis_BiomedicalSciences_MarckxNathan/assets/119006891/7553eb76-0b04-4767-bf4f-028d807a217d)
