# Mass Spectrometry Data Cleaning for Improved Identification of Post-Translational Modifications
## Data
### LC-MS data
For the development of the calibration module, .raw files from the following PRIDE projects was used and converted to mzML format:
- PXD001468 (https://www.ebi.ac.uk/pride/archive/projects/PXD001468)
- PXD032235 (https://www.ebi.ac.uk/pride/archive/projects/PXD032235)

The Python scripts to calibrate all reprocessed mzML files can be found in the "scripts" folder ("LRcal.py" and "RFcal.py").

### Search parameters
The following Sage search parameter JSON files were used and can be found in the "data" folder:
- Initial closed searches: SearchClosed.json
- Open search window selection: SearchClosed5ppm.json to SearchClosed30ppm.json
- Open searches: SearchOpen.json

### Sequence database
One FASTA file was consistently used through all the searches and can also be found in the "data" folder:
- human_ref_prot_and_contaminants.fasta

The custom Python script to combine the human sequences and contaminant sequences ("merge_fastas.py") can be found in the "scripts" folder.
The original human ("human_reference_proteome_2024_03_06.fasta") and contaminant ("GPM_crap.fasta") sequence FASTA files can be found in the "data" folder.

## Mass calibration pipeline
All Python an Jupyter Notebook files written for the construction of the mass calibration module can be found in the "scripts" folder.

The backbonde of the calibration pipeline itself was based on the 'Mass Spectrometry Workflow' Nextflow tool from GitHub user jvdnheme (https://github.ugent.be/jvdnheme/automatic_annotation). Extra modules including functionality to convert .xz into .raw files and to calibrate mzML files (LR and RF) were incorporated. The actual Nextflow script ("CalibScript.nf"), as well as a corresponding Nextflow script configuration file ("nextflow.config") can be found in the "CalibrationModule" folder

## Pipeline Schematic Representation
![pipeline](https://github.com/NathanMRX2404/Thesis_BiomedicalSciences_MarckxNathan/assets/119006891/7553eb76-0b04-4767-bf4f-028d807a217d)
