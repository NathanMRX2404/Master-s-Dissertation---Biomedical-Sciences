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

# RF models
import pandas as pd
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_squared_error

# Load the dataset
data = df_uncal_results_nonlin_q

# Get unique filenames
unique_filenames = np.sort(df_uncal_results_nonlin_q['filename'].unique())

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
for filename in np.sort(df_uncal_results_nonlin_q['filename'].unique()):
    # Filter dataframe for the current filename
    df_subset = df_uncal_results_nonlin_q[df_uncal_results_nonlin_q['filename'] == filename]
    
    # Get features for prediction
    X_subset = df_subset[['expMZ', 'RT', 'TIC', 'IT']]
    
    # Make predictions using the corresponding trained model
    rf_model = trained_models[filename]
    y_pred = rf_model.predict(X_subset)
    
    # Subtract RF 'predicted' deltaMZ from deltaMZ to compute RF adjusted deltaMZ
    df_uncal_results_nonlin_q.loc[df_uncal_results_nonlin_q['filename'] == filename, 'RF adj deltaMZ'] = \
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
    
    input_file = os.path.join("results", "TRFP", "PXD032235", f"{filename}")
    output_file = os.path.join("results", "TRFP", "PXD032235", f"{filename}_RFcal.mzML")
    
    # Apply the RF calibration to the mzML file
    with open(input_file, 'rb') as in_stream, open(output_file, 'wb') as out_stream:
        MzMLTransformer(in_stream, out_stream, lambda spectrum: RFcal(spectrum, rf_model)).write()
