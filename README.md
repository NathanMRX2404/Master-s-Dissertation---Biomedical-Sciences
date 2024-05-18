# Mass Spectrometry Data Cleaning for Improved Identification of Post-Translational Modifications
## Data
### LC-MS data
For the development of the calibration module, .raw files from the following PRIDE projects was used and converted to mzML format:
- PXD001468 (https://www.ebi.ac.uk/pride/archive/projects/PXD001468)
- PXD032235 (https://www.ebi.ac.uk/pride/archive/projects/PXD032235)
### Search parameters
The following Sage search parameter JSON files were used and can be found in the "data" folder:
Initial closed searches:
- SearchClosed.json
Closed searches for selecting open search fragment mass tolerance window:
- SearchClosed5ppm.json
- SearchClosed10ppm.json
- SearchClosed15ppm.json
- SearchClosed20ppm.json
- SearchClosed25ppm.json
- SearchClosed30ppm.json
Open searches:
- SearchOpen.json

## Mass Calibration Pipeline
The calibration pipeline was based on the Mass Spectrometry Workflow tool from GitHub user jvdnheme (https://github.ugent.be/jvdnheme/automatic_annotation)
The

## Pipeline Requirements and Installation
...
## Pipeline Schematic Representation
![pipeline](https://github.com/NathanMRX2404/Thesis_BiomedicalSciences_MarckxNathan/assets/119006891/7553eb76-0b04-4767-bf4f-028d807a217d)
