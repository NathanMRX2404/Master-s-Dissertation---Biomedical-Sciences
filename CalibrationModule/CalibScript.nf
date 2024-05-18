#!/usr/bin/env nextflow
//nextflow.enable.dsl=2

/*
 * Define the (standard) parameters
 */

// path to output .raw files after Decompress
params.output_path_decompress = "./data/PXD032235_Full"

// path to output of all processes after Sage
params.output_path = "./results"

// path to output of all processes after LRcalibration
params.output_path_LRcalibration = "./data"

// Decompress .raw.xz files
params.compressed_file_location = "/public/conode55_pride/PRIDE_DATA/PXD032235/RAW/*_HM_HEK_*"
params.run_decompress = false

// TRFP
params.raw_file_location = "data/PXD032235_Full/*.raw"
params.run_trfp = false

// TDF2MZML
params.tdf_folders_location = "data/*.d"
params.run_tdf2mzml = false

// SAGE
params.fasta_files_location = "data/human_ref_prot_and_contaminants.fasta"
params.sage_config_files_location = "data/sage/SearchOpen.json"
params.mzml_files_location = "results/TRFP/PXD032235/*_RFcal.mzML"
params.sage_max_forks = 2
params.splitFasta_size = 5.mb
params.run_sage = true

// LR calibration
params.run_LRcalibration = false

// ms2rescore
params.ms2rescore_config_files_location = "data/ms2rescore/*.json"
params.run_ms2rescore = false

// FlashLFQ
params.run_FlashLFQ = false

// organism_count_mass_tolerance
params.run_organism_count_mass_tolerance = false

/*
* Define the processes
*/

process decompressFiles {

    publishDir "${params.output_path_decompress}", mode: 'copy'

    input:
    path compressed_file 

    output:
    file "*.raw"

    script:
    """
    unxz --keep $compressed_file
    """
}

process runThermoRawFileParser{
    
    publishDir "${params.output_path}/TRFP", mode: 'copy'

    input:
    path raw_file 

    output:
    file "*.mzML"

    script:
    """
    thermorawfileparser --input=$raw_file -o=./ --format=1
    """

}

process runTDF2MZML{

    input:
    path d_file

    output:
    path "result.mzML"

    script:
    """
    tdf2mzml.py -i $d_file -o result.mzML
    """

}

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



process LRcalibration {

    publishDir "${params.output_path_LRcalibration}", mode: 'copy'

    input:
    path combined_sage_tsv_file

    output:
    path "${configFile}_${fastaFile}_SageResults"

    script:
    """
#!/usr/bin/env python3
import pandas as pd

# Read the combined Sage TSV file into a DataFrame
df = pd.read_csv("$combined_sage_tsv_file", sep="\t")

df_q = df[df["spectrum_q"] < 0.01]

# Save the filtered DataFrame to a new TSV file
df_q.to_csv("LRcalibrated_results.tsv", sep="\t", index=False)
    """
}


process ms2rescore {

    publishDir(
        path: "${params.output_path}/ms2rescore",
        mode: 'copy',
    )

    input:
    path configFile
    path sageResults
    path mzmlFiles
    path fastaFiles

    //output:
    //path "*.flashlfq.txt"

    script:
    """
    #!/bin/bash
    IFS=' ' read -r -a sageResults <<< "$sageResults"
    IFS=' ' read -r -a fastaFiles <<< "$fastaFiles"

    for sageResult in "\${sageResults[@]}"
    do
        fasta_name=\$(echo \$sageResult | cut -d'.' -f1)
        config_name=\$(echo \$sageResult | cut -d'.' -f2)
        
        for fastaFile in "\${fastaFiles[@]}"
        do
            if [ \${fastaFile%.*} == \$fasta_name ]; then
                ms2rescore -c $configFile -p \$sageResult -s . -f \$fastaFile -o "\${fasta_name}.\${config_name}.ms2rescore"
            fi
        done

    done
    """

    
}

process combine_sageTSV {

    publishDir(
        path: "${params.output_path}/sage",
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





process organism_count_mass_tolerance {
    
    publishDir(
        path: "${params.output_path}/organism_count_mass_tolerance",
        mode: 'copy',
    )


    input:
    path tsvFiles

    output:
    path "*_protein_count_organism.csv"
    path "*_ppm_distribution.png"

    script:
    """
    #!/usr/bin/env python3.8
    import pandas as pd
    from collections import Counter
    from matplotlib import pyplot as plt
    import numpy as np
    import sys

    _tsvFiles = "$tsvFiles"
    tsvFiles_list = _tsvFiles.split()

    for file in tsvFiles_list:

        df = pd.read_csv(file,sep="\t")

        file_name = file.replace(".tsv","")

        if len(df.index) == 0:
            sys.exit(-1)

        #filtered_df = df[df["peptide_q"] < 0.01]
        #filtered_df = filtered_df[filtered_df["label"] == 1]
        filtered_df = df[df["label"] == 1]

        if len(filtered_df.index) == 0:
            sys.exit(-2)

        filtered_df["organism"] = filtered_df["proteins"].str.split("_").str[1].str.split(";").str[0]

        count_organism = Counter(filtered_df["organism"])
        with open(file_name + "_protein_count_organism.csv", "w") as f:
            f.write("organism,count\\n")
            for k,v in  count_organism.most_common():
                f.write("{},{}\\n".format(k,v))
                

        mass_errors = filtered_df["expmass"]-filtered_df["calcmass"]

        ppm_errors = (mass_errors/filtered_df["expmass"])*1e6
        ppm_errors_filtered = ppm_errors[abs(mass_errors) < 0.5]
        filtered_mass_df = filtered_df[abs(mass_errors) < 0.5]

        plt.hist(ppm_errors_filtered,bins=500)
        try:
            plt.xlim(np.percentile(ppm_errors_filtered,90)*-1,np.percentile(ppm_errors_filtered,90))
        except IndexError:
            pass
        plt.savefig(file_name + "_ppm_distribution.png")
        plt.close()
    """
}

process runFlashLFQ {
    
    publishDir(
        path: "${params.output_path}/FlashLFQ ",
        mode: 'copy',
    )
    
    input:
    path flashLFQFile
    path mzmlFiles

    script:
    """
    #!/bin/bash
    echo task path: $PWD
    mkdir ${PWD}/mzmlFolder
    mv $mzmlFiles ${PWD}/mzmlFolder
    echo CMD.exe --idt "$flashLFQFile" --rep "${PWD}/mzmlFolder"
    """

}

/*
* Define the workflow
*/

workflow {

    /*
    * decompressFiles
    */

    def Decompress_output_ch

    if (params.run_decompress) {

        // Provide paths to the input raw.xz files for decompressFiles
        def input_compressed_files_ch = Channel.fromPath(params.compressed_file_location)

        // Execute the process decompressFiles with the input files and store the output in Decompress_output_ch channel
        Decompress_output_ch = decompressFiles(input_compressed_files_ch.flatten())
    }


    /*
    * runThermoRawFileParser
    */

    def TRFP_output_ch

    if (params.run_trfp) {
        if (params.run_decompress) {
            // Use the output of decompression as input for TRFP
            input_raw_files_ch = Decompress_output_ch
        } else {
            // Provide paths to the input raw files for TRFP
            input_raw_files_ch = Channel.fromPath(params.raw_file_location) 
        }
        // Execute the process runThermoRawFileParser with the input files and store the output in TRFP_output_ch channel
        TRFP_output_ch = runThermoRawFileParser(input_raw_files_ch.flatten())
    }

    /*
    * runTDF2MZML
    */

    def runTDF2MZML_output_ch

    if (params.run_tdf2mzml) {

        // Provide paths to the input files for runTDF2MZML
        def input_tdf_folders_ch = Channel.fromPath(params.tdf_folders_location, type: 'dir')

        // Execute the process runTDF2MZML with the input files and store the output in runTDF2MZML_output_ch channel
        runTDF2MZML_output_ch = runTDF2MZML(input_tdf_folders_ch)

    }
     
    /*
    * runSage
    */

    if (params.run_sage) {
        // Provide paths to the input files for runSage
        def input_fasta_ch = Channel.fromPath(params.fasta_files_location)
                                    
        def input_runSage_config_ch = Channel.fromPath(params.sage_config_files_location) 

        // mix the mzml files, output of runTDF2MZML and TRFP files
        def sage_input_mzml_ch = Channel.fromPath(params.mzml_files_location)

        if (params.run_trfp) {

            sage_input_mzml_ch = sage_input_mzml_ch.mix(TRFP_output_ch)

        }

        if (params.run_tdf2mzml) {

            sage_input_mzml_ch = sage_input_mzml_ch.mix(runTDF2MZML_output_ch)

        }

    // Execute the process runsage with input files and sage_input_mzml_ch
        def SAGE_output_ch = runSage(input_fasta_ch.splitFasta( size: params.splitFasta_size, file: true) ,input_runSage_config_ch,sage_input_mzml_ch.collect())

        def combine_sageTSV_output_ch = combine_sageTSV(SAGE_output_ch.collect())

    //combine_sageTSV_output_ch.view{"combine_sageTSV_output_ch: $it"}
    }



    /*
    * LRcalibration
    */

    if (params.run_LRcalibration) {

        // Provide paths to the input files for LRcalibration
        def input_LRcalibration_ch = combine_sageTSV_output_ch

        def LRcalibration_output_ch = LRcalibration(input_LRcalibration_ch)

    }   


    /*
    * ms2rescoregit
    */
    
    if (params.run_ms2rescore) {

        // Provide paths to the input files for ms2rescore
        def input_ms2rescore_config_ch = Channel.fromPath(params.ms2rescore_config_files_location)
        
        def ms2rescore_output_ch = ms2rescore(input_ms2rescore_config_ch, combine_sageTSV_output_ch, sage_input_mzml_ch.collect(),input_fasta_ch.collect())


        if (params.run_FlashLFQ) {
            runFlashLFQ(ms2rescore_output_ch, Channel.fromPath("data"))
        } 
    }

    /*
    * organism_count_mass_tolerance
    */

    if (params.run_organism_count_mass_tolerance) {
        organism_count_mass_tolerance(combine_sageTSV_output_ch)
    }

}
