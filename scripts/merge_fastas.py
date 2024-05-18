from Bio import SeqIO

def merge_fasta_files(input_file1, input_file2, output_file):
    records = []

    # Read the first input .fasta file
    for record in SeqIO.parse(input_file1, 'fasta'):
        records.append(record)

    # Read the second input .fasta file
    for record in SeqIO.parse(input_file2, 'fasta'):
        records.append(record)

    # Write the records to the output file
    SeqIO.write(records, output_file, 'fasta')

# Call the function with your input and output files
merge_fasta_files('human_reference_proteome_2024_03_06.fasta', 'GPM_crap.fasta', 'human_ref_prot_and_contaminants.fasta')
