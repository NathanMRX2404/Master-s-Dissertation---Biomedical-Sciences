from psims.transform.mzml import MzMLTransformer
import pandas as pd
import os
import pickle
from pyteomics import mass
import numpy as np

# Load datasets
PXD001468 = "results/sage/FULL_PXD001468.tsv"
df_uncal_results = pd.read_csv(PXD001468, sep='\t')
print("Lineair data uncalibrated CS:                     ", df_uncal_results.shape)
df_uncal_results_q = df_uncal_results[df_uncal_results["spectrum_q"] < 0.01]
print("Lineair data uncalibrated CS q:                   ", df_uncal_results_q.shape)

PXD032235 = "results/sage/FULL_PXD032235.tsv"
df_uncal_results_nonlin = pd.read_csv(PXD032235, sep='\t')
print("Lineair data uncalibrated CS:                     ", df_uncal_results_nonlin.shape)
df_uncal_results_nonlin_q = df_uncal_results_nonlin[df_uncal_results_nonlin["spectrum_q"] < 0.01]
print("Lineair data uncalibrated CS q:                   ", df_uncal_results_nonlin_q.shape)

df_uncal_results_q = df_uncal_results_q[['filename', 'scannr', 'peptide']]
df_uncal_results_nonlin_q = df_uncal_results_nonlin_q[['filename', 'scannr', 'peptide']]

# Load the dictionaries
with open("spectra_by_file.pkl", "rb") as f:
    spectra_by_file = pickle.load(f)

with open("spectra_by_file_nonlin.pkl", "rb") as f:
    spectra_by_file_nonlin = pickle.load(f)

# Initialize dictionary for mass calculations
modification_dict = {
    '[+15.9949]': 'ox',
    '[+57.0214]': 'cm'
}
print(modification_dict)

db = mass.Unimod()
aa_comp = dict(mass.std_aa_comp)
aa_comp['ox'] = db.by_title('Oxidation')['composition']
aa_comp['cm'] = db.by_title('Carbamidomethyl')['composition']
print(aa_comp)

# Construct 2 functions to extract features from linear and non-linear data
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

def extract_features_uncal_nonlin(row):
    file = row['filename']
    
    scannr_value = int(row['scannr'].split("scan=")[-1])
    index_in_mzml = scannr_value - 1
    
    spectrum = spectra_by_file_nonlin[file][index_in_mzml]
    
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

# Apply the functions to the dataframes
df_uncal_results_q[['expMZ', 'calcMZ', 'deltaMZ', 'RT', 'TIC', 'IT']] = df_uncal_results_q.apply(extract_features_uncal, axis=1, result_type='expand')
df_uncal_features = df_uncal_results_q[['expMZ', 'RT', 'TIC', 'IT', 'deltaMZ']]

df_uncal_results_nonlin_q[['expMZ', 'calcMZ', 'deltaMZ', 'RT', 'TIC', 'IT']] = df_uncal_results_nonlin_q.apply(extract_features_uncal_nonlin, axis=1, result_type='expand')
df_uncal_features_nonlin = df_uncal_results_nonlin_q[['expMZ', 'RT', 'TIC', 'IT', 'deltaMZ']]

# Construct datafram with filename, slope and intercept
# Get unique filenames and sort alphabetically
unique_filenames = np.sort(df_uncal_results_nonlin_q['filename'].unique())

# Loop through unique filenames and update the dataframe
for filename in unique_filenames:
    df_subset = df_uncal_results_nonlin_q[df_uncal_results_nonlin_q['filename'] == filename]
    
    x = df_subset['RT']
    y = df_subset['deltaMZ']
    coefficients = np.polyfit(x, y, 1)
    
    # Print coefficients
    slope, intercept = coefficients
    
    # Add coefficients to dataframe with increased precision
    df_uncal_results_nonlin_q.loc[df_uncal_results_nonlin_q['filename'] == filename, 'slope'] = slope
    df_uncal_results_nonlin_q.loc[df_uncal_results_nonlin_q['filename'] == filename, 'intercept'] = intercept
    
    # Subtract LR 'predicted' deltaMZ from deltaMZ to compute LR adjusted deltaMZ
    df_uncal_results_nonlin_q.loc[df_uncal_results_nonlin_q['filename'] == filename, 'LR adj deltaMZ'] = \
        df_subset['deltaMZ'] - (df_subset['RT'] * slope + intercept)

# Select the desired columns and keep the first occurrence of each unique filename
df_subset_LR_features = df_uncal_results_nonlin_q[['filename', 'slope', 'intercept']].groupby('filename').first().reset_index()
df_subset_LR_features

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
    input_file = os.path.join("results", "TRFP", "PXD032235", f"{filename}")
    output_file = os.path.join("results", "TRFP", "PXD032235", f"{filename}_LRcal.mzML")
    
    # Apply the LR calibration to the mzML file
    with open(input_file, 'rb') as in_stream, open(output_file, 'wb') as out_stream:
        MzMLTransformer(in_stream, out_stream, lambda spectrum: LRcal(spectrum, slope, intercept)).write()
