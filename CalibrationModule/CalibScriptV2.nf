#!/usr/bin/env nextflow
//nextflow.enable.dsl=2

// This is a Nextflow script originally adapted from GitHub user 'jvdnheme'. The file was modified to include 
// the necessary Python mass calibration pipeline.

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/*
 * Define the (standard) parameters
 */

// A) Raw file decompression
params.input_path_decompress = "/public/conode53_pride/PRIDE_DATA/PXD001468/RAW/b1948_293T_proteinID_12B_QE3_122212.raw.xz"
params.output_path_decompress = "./results/decompress"
params.run_decompress = false

// B1) TRFP
params.input_path_TRFP = "./results/decompress/20210115_HM_Ecoli_IAA_300minGr_R1.raw"
params.output_path_TRFP = "./results/TRFP"
params.run_trfp = false

// B2) TDF2MZML
params.input_path_TDF2MZML = "data/*.d"
params.output_path_TDF2MZML = "./results/TDF2MZML"
params.run_tdf2mzml = false

// C) Sage run
params.input_mzML_files_location = "./results/TRFP/*.mzML"
params.input_sage_config_files_location = "./data/sage/SearchClosed.json"
params.input_fasta_files_location = "./data/sage/human_ref_prot_and_contaminants.fasta"
params.output_path_sage = "./results/sage"
params.sage_max_forks = 2
params.splitFasta_size = 5.mb
params.run_sage = true

// D) Mass Calibration
params.output_path_mass_calibration = "./results/mass_calibration"
params.run_LR_mass_calibration = true
params.run_RF_mass_calibration = true

// E) Calibrated Sage run
params.run_calibrated_sage = true

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/*
* Define the processes
*/

// A) Decompress raw data (files from e.g. /public/conode folder can be directly decompressed)
process decompressFiles {
    
    if (params.run_decompress && !params.run_trfp && !params.run_tdf2mzml && !params.run_sage) {
        publishDir "${params.output_path_decompress}", mode: 'copy'
    }
    
    input:
    path compressed_file 

    output:
    file "*.raw"

    script:
    """
    unxz --keep $compressed_file
    """
}

// B1) ThermoRawFileParser
process runThermoRawFileParser {
    
    if ((!params.run_decompress && params.run_trfp && !params.run_tdf2mzml &&!params.run_sage) || (params.run_decompress && params.run_trfp && !params.run_tdf2mzml &&!params.run_sage)) {
        publishDir "${params.output_path_TRFP}", mode: 'copy'
    }    

    input:
    path raw_file 

    output:
    file "*.mzML"

    script:
    """
    thermorawfileparser --input=$raw_file -o=./ --format=1
    """
}

// B2) Topograph Data Files (TDF) to mzML
process runTDF2MZML{

    if (!params.run_decompress && !params.run_trfp && params.run_tdf2mzml &&!params.run_sage) {
        publishDir "${params.output_path_TDF2MZML}", mode: 'copy'
    } 

    input:
    path d_file

    output:
    path "result.mzML"

    script:
    """
    tdf2mzml.py -i $d_file -o result.mzML
    """
}

// C) Sage
process runSage {

    maxForks params.sage_max_forks
    errorStrategy 'retry'

    input:
    each path(fastaFile)
    each path(configFile)
    path mzmlFiles

    output:
    path "${configFile}_${fastaFile}_SageResults"

    script:
    """
    sage -f $fastaFile $configFile $mzmlFiles -o "${configFile}_${fastaFile}_SageResults"
    """
}

// C) Combine Sage TSV files
process combine_sageTSV {

    publishDir(
        path: "${params.output_path_sage}",
        mode: 'copy',
    )

    input:
    path sageFolders

    output:
    path "*.tsv"

    script:
    """
#!/usr/bin/env python3
import pandas as pd
import os

sageFolders = "$sageFolders"
sageFolders_list = sageFolders.split()

folders_map = {}

for folder_name in sageFolders_list:

    parts = folder_name.split(".")

    config_file = parts[0]
    fasta_file = parts[1]
    fasta_file = fasta_file.replace("json_", "")

    combination = fasta_file + "." + config_file
    if (combination) in folders_map:
        folders_map[combination].append(folder_name)
    else:
        folders_map[combination] = [folder_name]

for combination in folders_map:
    #print("combo:")
    #print(combination)
    dataframes = []
    for folder in folders_map[combination]:
        csv_file = [f for f in os.listdir(folder) if f.endswith(".tsv")][0]
        #print("FOLDER: {x}\\n FILE: {y}".format(x=folder, y=csv_file))
        dataframes.append(pd.read_csv(folder + "/" + csv_file, sep="\t"))

    full_df = pd.concat(dataframes)
    full_df.sort_values("sage_discriminant_score", ascending=False, inplace=True)
    full_df.drop_duplicates(subset=["scannr","filename"], inplace=True,keep="first")
    full_df.to_csv(combination + ".results.sage.tsv", sep="\t", index=False)
    
    """
}

// D)  LR Mass Calibration
process runMassCalibrationLR {

    publishDir(
        path: "${params.output_path_mass_calibration}",
        mode: 'copy',
    )

    input:
    path tsvFile
    path mzmlFiles

    output:
    path "*.mzML"

    script:
"""
#!/usr/bin/env python3
from psims.transform.mzml import MzMLTransformer
import pandas as pd
import os
import sys
import pickle
from pyteomics import mass
import numpy as np
from pyteomics import mzml

import pandas as pd
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_squared_error

# Read the TSV file into a DataFrame
tsv_file_path = "$tsvFile"
df_uncal_results = pd.read_csv(tsv_file_path, sep='\t')
df_uncal_results_q = df_uncal_results[df_uncal_results["spectrum_q"] < 0.01]
df_uncal_results_q = df_uncal_results_q[['filename', 'scannr', 'peptide']]

# Read the mzML files and store spectra into a dictionary
spectra_by_file = {}
for mzml_file in "$mzmlFiles".split():
    # Initialize lists for spectra and counters for each file
    spectra_list = []
        
    # Open the mzML file for reading
    with mzml.read(mzml_file) as reader:
        # Iterate over each spectrum in the mzML file
        for spectrum in reader:
            # Append the spectrum to the list
            spectra_list.append(spectrum)

    # Store spectra list and counter in dictionary
    spectra_by_file[mzml_file] = spectra_list

# Initialize dictionary for mass calculations
modification_dict = {
    '[+15.9949]': 'ox',
    '[+57.0214]': 'cm'
}
db = mass.Unimod()
aa_comp = dict(mass.std_aa_comp)
aa_comp['ox'] = db.by_title('Oxidation')['composition']
aa_comp['cm'] = db.by_title('Carbamidomethyl')['composition']

def extract_features_uncal(row):
    file = row['filename']
    
    scannr_value = int(row['scannr'].split("scan=")[-1])
    index_in_mzml = scannr_value - 1
    
    spectrum = spectra_by_file[file][index_in_mzml]
    
    # Check if 'precursorList' is present in the spectrum dictionary
    if 'precursorList' in spectrum:
        charge = spectrum['precursorList']['precursor'][0]['selectedIonList']['selectedIon'][0]['charge state']
    else:
        charge = 0
        
    sequence = row['peptide']

    #Get expMZ
    if 'precursorList' in spectrum and 'precursor' in spectrum['precursorList'] \
        and spectrum['precursorList']['precursor'] \
        and 'selectedIonList' in spectrum['precursorList']['precursor'][0] \
        and spectrum['precursorList']['precursor'][0]['selectedIonList']['selectedIon']:

        selected_ion_list = spectrum['precursorList']['precursor'][0]['selectedIonList']['selectedIon']
        for selected_ion in selected_ion_list:
            selected_ion_mz_str = selected_ion.get('selected ion m/z')
            if selected_ion_mz_str is not None:
                # Extract the numerical value from the string representation
                selected_ion_mz = float(selected_ion_mz_str)
                expMZ = float(selected_ion_mz)
                
                #return calcMZ, selected_ion_mz

    #get calcMZ
    for key, value in modification_dict.items():
        sequence = sequence.replace(key, value)
    seq = sequence
    if seq[-2:] == 'cm' or seq[-2:] == 'ox':
        if seq[-2:] == 'cm':
            extra_mass = mass.calculate_mass(aa_comp['cm'])
            seq = seq[:-2]
        elif seq[-2:] == 'ox':
            extra_mass = mass.calculate_mass(aa_comp['ox'])
            seq = seq[:-2]
    else:
        extra_mass = 0
    calmass = mass.calculate_mass(seq, aa_comp=aa_comp) + extra_mass
    if charge == 0 or charge is None:
        calcMZ = float(calmass)
        #return calcMZ
    else:
        calcMZ = float((calmass / charge) + 1.0072764667700085)
        #return calcMZ
        
    #Calculate deltaMZ
    deltaMZ = float(expMZ - calcMZ)
    # return expMZ, calcMZ, deltaMZ

    #get RT
    RT = float(spectrum['scanList']['scan'][0]['scan start time'])
    
    #get TIC
    TIC = float(spectrum['total ion current'])
        
    #get IT
    IT = float(spectrum['scanList']['scan'][0]['ion injection time'])
    
    return expMZ, calcMZ, deltaMZ, RT, TIC, IT

# Apply the function to the dataframes
df_uncal_results_q[['expMZ', 'calcMZ', 'deltaMZ', 'RT', 'TIC', 'IT']] = df_uncal_results_q.apply(extract_features_uncal, axis=1, result_type='expand')
#df_uncal_features = df_uncal_results_q[['expMZ', 'RT', 'TIC', 'IT', 'deltaMZ']]

# Get unique filenames and sort alphabetically
unique_filenames = np.sort(df_uncal_results_q['filename'].unique())

# Loop through unique filenames and update the dataframe
for filename in unique_filenames:
    df_subset = df_uncal_results_q[df_uncal_results_q['filename'] == filename]
    
    x = df_subset['RT']
    y = df_subset['deltaMZ']
    coefficients = np.polyfit(x, y, 1)
    
    # Print coefficients
    slope, intercept = coefficients
    
    # Add coefficients to dataframe with increased precision
    df_uncal_results_q.loc[df_uncal_results_q['filename'] == filename, 'slope'] = slope
    df_uncal_results_q.loc[df_uncal_results_q['filename'] == filename, 'intercept'] = intercept
    
    # Subtract LR 'predicted' deltaMZ from deltaMZ to compute LR adjusted deltaMZ
    df_uncal_results_q.loc[df_uncal_results_q['filename'] == filename, 'LR adj deltaMZ'] = \
        df_subset['deltaMZ'] - (df_subset['RT'] * slope + intercept)

# ----- LR calibration
# Select the desired columns and keep the first occurrence of each unique filename
df_subset_LR_features = df_uncal_results_q[['filename', 'slope', 'intercept']].groupby('filename').first().reset_index()

# Function for LR calibration
def LRcal(spectrum, slope, intercept):
    retention_time = spectrum['scanList']['scan'][0]['scan start time']
    mz_values = spectrum['m/z array']
    
    # Apply the calibration function to the m/z values
    calibrated_mz_values = mz_values - (retention_time * slope + intercept)
    
    # Replace the original m/z values with the calibrated ones
    spectrum['m/z array'] = calibrated_mz_values   
    
    # If there are selected ions, calibrate their m/z values as well
    if 'precursorList' in spectrum and 'precursor' in spectrum['precursorList'] \
            and spectrum['precursorList']['precursor'] \
            and 'selectedIonList' in spectrum['precursorList']['precursor'][0] \
            and spectrum['precursorList']['precursor'][0]['selectedIonList']['selectedIon']:
        
        selected_ion_list = spectrum['precursorList']['precursor'][0]['selectedIonList']['selectedIon']
        for selected_ion in selected_ion_list:
            selected_ion_mz = selected_ion.get('selected ion m/z')
            if selected_ion_mz is not None:
                selected_ion_mz_calibrated = selected_ion_mz - (retention_time * slope + intercept)
                selected_ion['selected ion m/z'] = selected_ion_mz_calibrated
    
    return spectrum


# Loop through each mzML file and apply calibration
for filename in df_subset_LR_features['filename']:
    # Get the slope and intercept for the current filename
    slope = df_subset_LR_features[df_subset_LR_features['filename'] == filename]['slope'].values[0]
    intercept = df_subset_LR_features[df_subset_LR_features['filename'] == filename]['intercept'].values[0]
    
    # Define input and output file paths


    input_file = filename
    output_file = filename.replace(".mzML", "_LRcal.mzML")
    
    # Apply the LR calibration to the mzML file
    with open(input_file, 'rb') as in_stream, open(output_file, 'wb') as out_stream:
        MzMLTransformer(in_stream, out_stream, lambda spectrum: LRcal(spectrum, slope, intercept)).write()

"""
}

// D)  RF Mass Calibration
process runMassCalibrationRF {

    publishDir(
        path: "${params.output_path_mass_calibration}",
        mode: 'copy',
    )

    input:
    path tsvFile
    path mzmlFiles

    output:
    path "*.mzML"

    script:
"""
#!/usr/bin/env python3
from psims.transform.mzml import MzMLTransformer
import pandas as pd
import os
import sys
import pickle
from pyteomics import mass
import numpy as np
from pyteomics import mzml

import pandas as pd
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_squared_error

# Read the TSV file into a DataFrame
tsv_file_path = "$tsvFile"
df_uncal_results = pd.read_csv(tsv_file_path, sep='\t')
df_uncal_results_q = df_uncal_results[df_uncal_results["spectrum_q"] < 0.01]
df_uncal_results_q = df_uncal_results_q[['filename', 'scannr', 'peptide']]

# Read the mzML files and store spectra into a dictionary
spectra_by_file = {}
for mzml_file in "$mzmlFiles".split():
    # Initialize lists for spectra and counters for each file
    spectra_list = []
        
    # Open the mzML file for reading
    with mzml.read(mzml_file) as reader:
        # Iterate over each spectrum in the mzML file
        for spectrum in reader:
            # Append the spectrum to the list
            spectra_list.append(spectrum)

    # Store spectra list and counter in dictionary
    spectra_by_file[mzml_file] = spectra_list

# Initialize dictionary for mass calculations
modification_dict = {
    '[+15.9949]': 'ox',
    '[+57.0214]': 'cm'
}
db = mass.Unimod()
aa_comp = dict(mass.std_aa_comp)
aa_comp['ox'] = db.by_title('Oxidation')['composition']
aa_comp['cm'] = db.by_title('Carbamidomethyl')['composition']

def extract_features_uncal(row):
    file = row['filename']
    
    scannr_value = int(row['scannr'].split("scan=")[-1])
    index_in_mzml = scannr_value - 1
    
    spectrum = spectra_by_file[file][index_in_mzml]
    
    # Check if 'precursorList' is present in the spectrum dictionary
    if 'precursorList' in spectrum:
        charge = spectrum['precursorList']['precursor'][0]['selectedIonList']['selectedIon'][0]['charge state']
    else:
        charge = 0
        
    sequence = row['peptide']

    #Get expMZ
    if 'precursorList' in spectrum and 'precursor' in spectrum['precursorList'] \
        and spectrum['precursorList']['precursor'] \
        and 'selectedIonList' in spectrum['precursorList']['precursor'][0] \
        and spectrum['precursorList']['precursor'][0]['selectedIonList']['selectedIon']:

        selected_ion_list = spectrum['precursorList']['precursor'][0]['selectedIonList']['selectedIon']
        for selected_ion in selected_ion_list:
            selected_ion_mz_str = selected_ion.get('selected ion m/z')
            if selected_ion_mz_str is not None:
                # Extract the numerical value from the string representation
                selected_ion_mz = float(selected_ion_mz_str)
                expMZ = float(selected_ion_mz)
                
                #return calcMZ, selected_ion_mz

    #get calcMZ
    for key, value in modification_dict.items():
        sequence = sequence.replace(key, value)
    seq = sequence
    if seq[-2:] == 'cm' or seq[-2:] == 'ox':
        if seq[-2:] == 'cm':
            extra_mass = mass.calculate_mass(aa_comp['cm'])
            seq = seq[:-2]
        elif seq[-2:] == 'ox':
            extra_mass = mass.calculate_mass(aa_comp['ox'])
            seq = seq[:-2]
    else:
        extra_mass = 0
    calmass = mass.calculate_mass(seq, aa_comp=aa_comp) + extra_mass
    if charge == 0 or charge is None:
        calcMZ = float(calmass)
        #return calcMZ
    else:
        calcMZ = float((calmass / charge) + 1.0072764667700085)
        #return calcMZ
        
    #Calculate deltaMZ
    deltaMZ = float(expMZ - calcMZ)
    # return expMZ, calcMZ, deltaMZ

    #get RT
    RT = float(spectrum['scanList']['scan'][0]['scan start time'])
    
    #get TIC
    TIC = float(spectrum['total ion current'])
        
    #get IT
    IT = float(spectrum['scanList']['scan'][0]['ion injection time'])
    
    return expMZ, calcMZ, deltaMZ, RT, TIC, IT

# Apply the function to the dataframes
df_uncal_results_q[['expMZ', 'calcMZ', 'deltaMZ', 'RT', 'TIC', 'IT']] = df_uncal_results_q.apply(extract_features_uncal, axis=1, result_type='expand')
#df_uncal_features = df_uncal_results_q[['expMZ', 'RT', 'TIC', 'IT', 'deltaMZ']]

# Get unique filenames and sort alphabetically
unique_filenames = np.sort(df_uncal_results_q['filename'].unique())

# Loop through unique filenames and update the dataframe
for filename in unique_filenames:
    df_subset = df_uncal_results_q[df_uncal_results_q['filename'] == filename]
    
    x = df_subset['RT']
    y = df_subset['deltaMZ']
    coefficients = np.polyfit(x, y, 1)
    
    # Print coefficients
    slope, intercept = coefficients
    
    # Add coefficients to dataframe with increased precision
    df_uncal_results_q.loc[df_uncal_results_q['filename'] == filename, 'slope'] = slope
    df_uncal_results_q.loc[df_uncal_results_q['filename'] == filename, 'intercept'] = intercept
    
    # Subtract LR 'predicted' deltaMZ from deltaMZ to compute LR adjusted deltaMZ
    df_uncal_results_q.loc[df_uncal_results_q['filename'] == filename, 'LR adj deltaMZ'] = \
        df_subset['deltaMZ'] - (df_subset['RT'] * slope + intercept)

# ----- RF calibration
# Load the dataset
data = df_uncal_results_q

# Get unique filenames
unique_filenames = np.sort(df_uncal_results_q['filename'].unique())

# Dictionary to store trained models
trained_models = {}

# Iterate over unique filenames
for filename in unique_filenames:
    # Filter data for the current filename
    subset_data = data[data['filename'] == filename]
    
    # Split the subset dataset into features (X) and labels (y)
    X_subset = subset_data[['expMZ', 'RT', 'TIC', 'IT']]  # Features
    y_subset = subset_data['deltaMZ']  # Labels
    
    # Split the dataset into training and testing sets
    X_train, X_test, y_train, y_test = train_test_split(X_subset, y_subset, test_size=0.2, random_state=42)
    
    # Train the Random Forest Regressor model
    rf_regressor = RandomForestRegressor(n_estimators=100, random_state=42)
    rf_regressor.fit(X_subset, y_subset)
    
    # Store trained model in the dictionary
    trained_models[filename] = rf_regressor
    
    # Optionally, print the MSE for each subset
    y_pred_subset = rf_regressor.predict(X_test)
    mse_subset = mean_squared_error(y_test, y_pred_subset)
    print(f"Filename: {filename}, Mean Squared Error (MSE): {mse_subset}")



# Loop through unique filenames to compute LR adjusted deltaMZ
for filename in np.sort(df_uncal_results_q['filename'].unique()):
    # Filter dataframe for the current filename
    df_subset = df_uncal_results_q[df_uncal_results_q['filename'] == filename]
    
    # Get features for prediction
    X_subset = df_subset[['expMZ', 'RT', 'TIC', 'IT']]
    
    # Make predictions using the corresponding trained model
    rf_model = trained_models[filename]
    y_pred = rf_model.predict(X_subset)
    
    # Subtract RF 'predicted' deltaMZ from deltaMZ to compute RF adjusted deltaMZ
    df_uncal_results_q.loc[df_uncal_results_q['filename'] == filename, 'RF adj deltaMZ'] = \
        df_subset['deltaMZ'] - y_pred

import warnings
import numpy as np
from psims.transform.mzml import MzMLTransformer

# Suppress UserWarning about feature names
warnings.filterwarnings("ignore", category=UserWarning)

def get_x(spectrum):
    mz_values = spectrum['m/z array']
    retention_time = float(spectrum['scanList']['scan'][0]['scan start time'])
    total_ion_count = float(spectrum['total ion current'])
    injection_time = float(spectrum['scanList']['scan'][0]['ion injection time'])
    x = [[mz_value, retention_time, total_ion_count, injection_time] for mz_value in mz_values]
    return x

def RFcal(spectrum, rf_model):
    mz_values = spectrum['m/z array']
    x = get_x(spectrum)
    retention_time = float(spectrum['scanList']['scan'][0]['scan start time'])
    total_ion_count = float(spectrum['total ion current'])
    injection_time = float(spectrum['scanList']['scan'][0]['ion injection time'])
    
    # Apply the calibration function to the m/z values
    calibrated_mz_values = np.array(mz_values) - np.array(rf_model.predict(x))
    
    # Replace the original m/z values with the calibrated ones
    spectrum['m/z array'] = calibrated_mz_values   
    
    # If there are selected ions, calibrate their m/z values as well
    if 'precursorList' in spectrum and 'precursor' in spectrum['precursorList'] \
            and spectrum['precursorList']['precursor'] \
            and 'selectedIonList' in spectrum['precursorList']['precursor'][0] \
            and spectrum['precursorList']['precursor'][0]['selectedIonList']['selectedIon']:
        
        selected_ion_list = spectrum['precursorList']['precursor'][0]['selectedIonList']['selectedIon']
        for selected_ion in selected_ion_list:
            selected_ion_mz_str = selected_ion.get('selected ion m/z')
            if selected_ion_mz_str is not None:
                # Extract the numerical value from the string representation
                selected_ion_mz = float(spectrum['precursorList']['precursor'][0]['selectedIonList']['selectedIon'][0]['selected ion m/z'])
                
                # Create features for prediction for the selected ion
                x_selected_ion = np.array([[selected_ion_mz, retention_time, total_ion_count, injection_time]])
                
                # Calibrate the selected ion m/z value
                calibrated_selected_ion_mz = str(selected_ion_mz - rf_model.predict(x_selected_ion))[1:-1]
                selected_ion['selected ion m/z'] = calibrated_selected_ion_mz
    
    return spectrum



import os
from psims.transform.mzml import MzMLTransformer

for filename in unique_filenames:
    
    rf_model = trained_models[filename]

    input_file = filename
    output_file = filename.replace(".mzML", "_RFcal.mzML")

    # Apply the RF calibration to the mzML file
    with open(input_file, 'rb') as in_stream, open(output_file, 'wb') as out_stream:
        MzMLTransformer(in_stream, out_stream, lambda spectrum: RFcal(spectrum, rf_model)).write()
"""
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/*
* Define the workflow
*/

workflow {

    /*
    * A) decompressFiles
    */

    def Decompress_output_ch

    if (params.run_decompress) {

        // Provide paths to the input raw.xz files for decompressFiles
        def input_compressed_files_ch = Channel.fromPath(params.input_path_decompress)

        // Execute the process decompressFiles with the input files and store the output in Decompress_output_ch channel
        Decompress_output_ch = decompressFiles(input_compressed_files_ch.flatten())
    }

    /*
    * B1) ThermoRawFileParser
    */

    def TRFP_output_ch

    if (params.run_trfp) {
        if (params.run_decompress) {
            // Use the output of decompression as input for TRFP
            input_raw_files_ch = Decompress_output_ch
        } else {
            // Provide paths to the input raw files for TRFP
            input_raw_files_ch = Channel.fromPath(params.input_path_TRFP) 
        }
        // Execute the process runThermoRawFileParser with the input files and store the output in TRFP_output_ch channel
        TRFP_output_ch = runThermoRawFileParser(input_raw_files_ch.flatten())
    }

    /*
    * B2) TDF2MZML
    */

    def runTDF2MZML_output_ch

    if (params.run_tdf2mzml) {

        // Provide paths to the input files for runTDF2MZML
        def input_tdf_folders_ch = Channel.fromPath(params.input_path_TDF2MZML, type: 'dir')

        // Execute the process runTDF2MZML with the input files and store the output in runTDF2MZML_output_ch channel
        runTDF2MZML_output_ch = runTDF2MZML(input_tdf_folders_ch)

    }

    /*
    * C) Sage
    */

    def combine_sageTSV_output_ch

    if (params.run_sage) {

        def sage_input_mzml_ch = Channel.empty()

        // Provide paths to the input files for runSage
        def input_fasta_ch = Channel.fromPath(params.input_fasta_files_location)
                                    
        def input_runSage_config_ch = Channel.fromPath(params.input_sage_config_files_location) 

        if (!params.run_trfp && !params.run_tdf2mzml) {
            // mix the mzml files, output of runTDF2MZML and TRFP files
            sage_input_mzml_ch = Channel.fromPath(params.input_mzML_files_location)
        }

        if (params.run_trfp || (params.run_decompress && params.run_trfp)) {

            sage_input_mzml_ch = sage_input_mzml_ch.mix(TRFP_output_ch)

        }

        if (params.run_tdf2mzml) {

            sage_input_mzml_ch = sage_input_mzml_ch.mix(runTDF2MZML_output_ch)

        }

    // Execute the process runsage with input files and sage_input_mzml_ch
        def SAGE_output_ch = runSage(input_fasta_ch.splitFasta( size: params.splitFasta_size, file: true) ,input_runSage_config_ch,sage_input_mzml_ch.collect())

        combine_sageTSV_output_ch = combine_sageTSV(SAGE_output_ch.collect())

    }

    /*
    * D) Mass Calibration
    */

    sage_input_mzml_ch = Channel.fromPath(params.input_mzML_files_location)

    if (params.run_LR_mass_calibration) {
            def runMassCalibrationLR_output_ch = runMassCalibrationLR(combine_sageTSV_output_ch, sage_input_mzml_ch)
    }

    if (params.run_RF_mass_calibration) {
            def runMassCalibrationRF_output_ch = runMassCalibrationRF(combine_sageTSV_output_ch, sage_input_mzml_ch)
    }

    /*
    * E) Calibrated Sage
    */

    if(params.run_calibrated_sage) {
        def input_fasta_ch = Channel.fromPath(params.input_fasta_files_location)
        def input_runSage_config_ch = Channel.fromPath(params.input_sage_config_files_location) 
        def calibrated_mzml_ch = Channel.empty()

        if (params.run_LR_mass_calibration && !params.run_RF_mass_calibration) {
            calibrated_mzml_ch = calibrated_mzml_ch.mix(runMassCalibrationLR_output_ch)
        }

        if (params.run_RF_mass_calibration && !params.run_LR_mass_calibration) {
            calibrated_mzml_ch = calibrated_mzml_ch.mix(runMassCalibrationRF_output_ch)
        }

        if (params.run_LR_mass_calibration && params.run_RF_mass_calibration) {
            calibrated_mzml_ch = calibrated_mzml_ch.mix(runMassCalibrationLR_output_ch, runMassCalibrationRF_output_ch)
        }

        def SAGE_Calibrated_output_ch = runSage(input_fasta_ch.splitFasta(size: params.splitFasta_size, file: true), input_runSage_config_ch, calibrated_mzml_ch.collect())
        combine_sageTSV(SAGE_Calibrated_output_ch.collect())
        }
    }
