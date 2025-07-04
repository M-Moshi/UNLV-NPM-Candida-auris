import pandas as pd
import sys

def extract_vcf_data(line):
    data = line.split('\t')
    if len(data) > 9:  # Check if FORMAT and sample data exist
        format_parts = data[8].split(':')
        sample_parts = data[9].split(':')
        format_dict = {f: sample_parts[i] for i, f in enumerate(format_parts) if i < len(sample_parts)}

        return [format_dict.get('REF_DP', ''), format_dict.get('ALT_DP', ''), format_dict.get('ALT_FREQ', '')]
    else:
        return ['', '', '']

def main(vcf_file_path, tsv_file_path, output_file_path):
    # Read the VCF file, skipping metadata lines
    with open(vcf_file_path, 'r') as file:
        vcf_data = [line.strip() for line in file if not line.startswith('#')]

    # Extract the data
    extracted_data = [extract_vcf_data(line) for line in vcf_data]

    # Create DataFrame from the VCF data
    vcf_df = pd.DataFrame(extracted_data, columns=['REF_DP', 'ALT_DP', 'ALT_FREQ'])

    # Read the TSV file generated by SnpSift
    tsv_df = pd.read_csv(tsv_file_path, sep='\t')

    # Combine the dataframes
    combined_df = pd.concat([tsv_df, vcf_df], axis=1)

    # Save the combined dataframe to a new TSV file
    combined_df.to_csv(output_file_path, sep='\t', index=False)

if __name__ == '__main__':
    if len(sys.argv) != 4:
        print("Usage: python process_vcf.py <vcf_file_path> <tsv_file_path> <output_file_path>")
        sys.exit(1)

    main(sys.argv[1], sys.argv[2], sys.argv[3])
