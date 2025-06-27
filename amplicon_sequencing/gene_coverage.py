import pandas as pd
import os
import glob
import argparse # Use argparse for robust command-line argument handling

def calculate_coverage(results_folder, positions_file):
    """
    Calculates the percentage of coverage for specific genomic regions based on depth files.

    Args:
        results_folder (str): The path to the main results directory. 
                              This directory must contain a 'depth/' subfolder.
        positions_file (str): The path to the TSV file containing the start and end positions.
    """
    
    # --- 1. Define Paths and Check for Inputs ---
    
    # Construct the path to the 'depth' folder automatically
    depth_folder = os.path.join(results_folder, 'depth')

    # Check if the required directories and files exist
    if not os.path.isdir(results_folder):
        print(f"Error: Results folder not found at '{results_folder}'")
        return
    if not os.path.isdir(depth_folder):
        print(f"Error: 'depth' subfolder not found inside '{results_folder}'")
        return
    if not os.path.isfile(positions_file):
        print(f"Error: Positions file not found at '{positions_file}'")
        return

    print(f"Results Folder: {results_folder}")
    print(f"Positions File: {positions_file}")

    # Create a new directory for the output reports inside the results folder
    output_dir = os.path.join(results_folder, 'coverage_analysis')
    os.makedirs(output_dir, exist_ok=True)
    print(f"Output will be saved to: {output_dir}")

    # --- 2. Load and Process Input Files ---

    # Load the positions TSV
    positions_df = pd.read_csv(positions_file, sep='\t')
    
    # Initialize a list to store the raw output data
    output_data = []

    # Get a list of all .udepth files in the depth folder
    udepth_files = glob.glob(os.path.join(depth_folder, '*.udepth'))
    print(f'\nFound {len(udepth_files)} udepth files to process.')

    if not udepth_files:
        print("Warning: No .udepth files found. Exiting.")
        return

    # --- 3. Main Processing Loop ---

    # Process each udepth file
    for udepth_file in udepth_files:
        sample_name = os.path.basename(udepth_file).replace('.udepth', '')
        print(f'--> Processing sample: {sample_name}')
        
        # Load the udepth file
        udepth_df = pd.read_csv(udepth_file, sep='\t', header=None, names=['Pos', 'Depth'])

        # Process each region defined in the positions file
        for idx, row in positions_df.iterrows():
            gene_name = row['File'] # Assuming 'File' column contains the gene/region name
            start_pos = row['Start']
            end_pos = row['End']
            
            # Filter the depth data for the current region
            filtered_udf = udepth_df[(udepth_df['Pos'] >= start_pos) & (udepth_df['Pos'] <= end_pos)]

            # Calculate the percentage of positions with depth > 50 and > 100
            if len(filtered_udf) > 0:
                percentage_over_50 = (filtered_udf['Depth'] > 50).mean() * 100
                percentage_over_100 = (filtered_udf['Depth'] > 100).mean() * 100
                output_data.append([gene_name, sample_name, percentage_over_50, percentage_over_100])
            else:
                # Handle cases where a region has no coverage data
                output_data.append([gene_name, sample_name, 0.0, 0.0])

    # --- 4. Generate and Save Output Reports ---

    # Convert the collected data to a DataFrame
    output_df = pd.DataFrame(output_data, columns=['Gene', 'Sample', 'Percent_Coverage_>50x', 'Percent_Coverage_>100x'])
    raw_output_path = os.path.join(output_dir, 'raw_coverage_report.tsv')
    output_df.to_csv(raw_output_path, sep='\t', index=False)
    print(f'\nRaw coverage report saved to: {raw_output_path}')

    # Create a pivot table for a more readable, matrix-style report
    # Using >50x coverage for the pivot table as an example
    pivot_df = output_df.pivot(index='Gene', columns='Sample', values='Percent_Coverage_>50x')
    pivot_output_path = os.path.join(output_dir, 'pivot_coverage_report_>50x.tsv')
    pivot_df.to_csv(pivot_output_path, sep='\t')
    print(f'Pivot table report saved to: {pivot_output_path}')

    print('\nAnalysis complete.')


if __name__ == '__main__':
    # --- Command-Line Argument Parsing ---
    parser = argparse.ArgumentParser(description="Calculate coverage percentage for genomic regions from sequencing depth files.")
    
    parser.add_argument("results_folder", 
                        type=str, 
                        help="Path to the main results folder containing the 'depth/' subdirectory.")
    
    parser.add_argument("positions_file", 
                        type=str, 
                        help="Path to the TSV file with 'File', 'Start', and 'End' columns for genomic regions.")
    
    args = parser.parse_args()
    
    # --- Run the main function ---
    calculate_coverage(args.results_folder, args.positions_file)
