import os
from pyteomics import mzml
import pickle

folder_path = "results/TRFP/PXD001468"

spectra_by_file_LR = {}

for file_name in os.listdir(folder_path):
    if file_name.endswith('_LRcal.mzML'):
        file_path = os.path.join(folder_path, file_name)
        
        # Initialize lists for spectra and counters for each file
        spectra_list = []
        
        # Open the mzML file for reading
        with mzml.read(file_path) as reader:
            # Iterate over each spectrum in the mzML file
            for spectrum in reader:
                # Append the spectrum to the list
                spectra_list.append(spectrum)
        
        # Store spectra list and counter in dictionary
        spectra_by_file_LR[file_name] = spectra_list
        
        print(f"File: {file_name}, Total spectra read: {len(spectra_list)}")
        
# Save the dictionary to a file
with open("spectra_by_file_LR.pkl", "wb") as f:
    pickle.dump(spectra_by_file_LR, f)
    
spectra_by_file_RF = {}

for file_name in os.listdir(folder_path):
    if file_name.endswith('_RFcal.mzML'):
        file_path = os.path.join(folder_path, file_name)
        
        # Initialize lists for spectra and counters for each file
        spectra_list = []
        
        # Open the mzML file for reading
        with mzml.read(file_path) as reader:
            # Iterate over each spectrum in the mzML file
            for spectrum in reader:
                # Append the spectrum to the list
                spectra_list.append(spectrum)
        
        # Store spectra list and counter in dictionary
        spectra_by_file_RF[file_name] = spectra_list
        
        print(f"File: {file_name}, Total spectra read: {len(spectra_list)}")
        
# Save the dictionary to a file
with open("spectra_by_file_RF.pkl", "wb") as f:
    pickle.dump(spectra_by_file_RF, f)
    
folder_path = "results/TRFP/PXD032235"   

spectra_by_file_LR_nonlin = {}

for file_name in os.listdir(folder_path):
    if file_name.endswith('_LRcal.mzML'):
        file_path = os.path.join(folder_path, file_name)
        
        # Initialize lists for spectra and counters for each file
        spectra_list = []
        
        # Open the mzML file for reading
        with mzml.read(file_path) as reader:
            # Iterate over each spectrum in the mzML file
            for spectrum in reader:
                # Append the spectrum to the list
                spectra_list.append(spectrum)
        
        # Store spectra list and counter in dictionary
        spectra_by_file_LR_nonlin[file_name] = spectra_list
        
        print(f"File: {file_name}, Total spectra read: {len(spectra_list)}")
        
# Save the dictionary to a file
with open("spectra_by_file_LR_nonlin.pkl", "wb") as f:
    pickle.dump(spectra_by_file_LR_nonlin, f)
    
spectra_by_file_RF_nonlin = {}

for file_name in os.listdir(folder_path):
    if file_name.endswith('_RFcal.mzML'):
        file_path = os.path.join(folder_path, file_name)
        
        # Initialize lists for spectra and counters for each file
        spectra_list = []
        
        # Open the mzML file for reading
        with mzml.read(file_path) as reader:
            # Iterate over each spectrum in the mzML file
            for spectrum in reader:
                # Append the spectrum to the list
                spectra_list.append(spectrum)
        
        # Store spectra list and counter in dictionary
        spectra_by_file_RF_nonlin[file_name] = spectra_list
        
        print(f"File: {file_name}, Total spectra read: {len(spectra_list)}")
        
# Save the dictionary to a file
with open("spectra_by_file_RF_nonlin.pkl", "wb") as f:
    pickle.dump(spectra_by_file_RF_nonlin, f)
