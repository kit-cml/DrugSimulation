#!/usr/bin/env python3

import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path
import re
import sys

# Helper function to check if a string can be converted to float
def _is_valid_float(s):
    try:
        float(s)
        return True
    except ValueError:
        return False

def extract_core_number(filename):
    """Extract core number from filename using regex"""
    match = re.search(r'core(\d+)', filename)
    return int(match.group(1)) if match else -1

def read_feature_data(base_path, drug_name, concentration):
    """
    Read all feature_core files for a given concentration
    Returns a sorted DataFrame with all features
    """
    feature_files = []
    folder_path = Path(base_path) / str(concentration)
    print(folder_path)
    
    # Find all feature files matching the pattern
    pattern = f"{drug_name}_{concentration}_features_core*.csv"
    for file in folder_path.glob(pattern):
        feature_files.append(str(file))
    
    if not feature_files:
        print(f"No feature files found for concentration {concentration}")
        return None
    
    # Read and combine all files
    dfs = []
    for file in sorted(feature_files, key=extract_core_number):
        try:
            df = pd.read_csv(file)
            df['concentration'] = concentration  # Add concentration as a column
            dfs.append(df)
        except Exception as e:
            print(f"Error reading {file}: {e}")
    
    if not dfs:
        return None
    
    # Combine all DataFrames
    combined_df = pd.concat(dfs, ignore_index=True)
    return combined_df

def create_feature_plots(all_data, drug_name, output_folder):
    """
    Create box plots for each feature showing all concentrations
    """
    # Make sure output folder exists
    os.makedirs(output_folder, exist_ok=True)
    
    # List of features to plot
    features = [
        'qNet', 'qInward', 'INaLAUC', 'ICaLAUC',
        'APD90', 'APD50', 'APDtri', 'Vmmax', 'Vmrest',
        'dVmdtmax', 'dVmdtmaxrepol', 'CaTD90', 'CaTD50',
        'CaTDtri', 'Camax', 'Carest'
    ]
    
    # Dictionary for feature name formatting
    feature_format = {
        'qNet': 'qNet (C/F)',
        'qInward': 'qInward (C/F)',
        'INaLAUC': 'INaL_AUC (A)',
        'ICaLAUC': 'ICaL_AUC (A)',
        'APD0': 'APD90 (ms)',
        'APD50': 'APD50 (ms)',
        'APDtri': 'APD_tri (ms)',
        'Vmmax': 'Vm_max (mV)',
        'Vmrest': 'Vm_rest (mV)',
        'dVmdtmax': 'dVm/dt_max (mV/ms)',
        'dVmdtmaxrepol': 'dVm/dt_max_repol (mV/ms)',
        'CaTD90': 'CaTD90 (ms)',
        'CaTD50': 'CaTD50 (ms)',
        'CaTDtri': 'CaTD_tri (ms)',
        'Capeak': 'Ca_peak (nM)',
        'Cavalley': 'Ca_valley (nM)'
    }
    
    # Get unique concentrations
    concentrations = sorted(all_data['concentration'].unique(), key=float)
    
    # Create a box plot for each feature
    for feature in features:
        fig, ax = plt.subplots(figsize=(8, 4))
        
        # Prepare data for boxplot
        data_to_plot = [all_data[all_data['concentration'] == conc][feature].values 
                       for conc in concentrations]
        
        # Create box plot
        box_plot = ax.boxplot(
            data_to_plot,
            patch_artist=False,
            medianprops=dict(color="black", linewidth=1.0),
            boxprops=dict(color='black', linewidth=1.0),
            whiskerprops=dict(color='black', linewidth=1.0),
            capprops=dict(color='black', linewidth=1.0),
            flierprops=dict(marker='o', markerfacecolor='black', markersize=4, 
                          markeredgecolor='black'),
            showfliers=True,
            whis=1.5
        )
        
        # Get formatted feature name
        formatted_feature = feature_format.get(feature, feature)
        
        # Remove top and right spines
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        
        # Set x-axis labels
        ax.set_xticks(range(1, len(concentrations) + 1))
        ax.set_xticklabels(concentrations)
        
        #plt.title(f'{formatted_feature}')
        plt.xlabel('Concentrations (nM)')
        plt.ylabel(formatted_feature)
        
        # Remove grid
        plt.grid(False)
        
        # Tight layout
        plt.tight_layout()
        
        # Save the plot
        plt.savefig(os.path.join(output_folder, f'{feature}_boxplot.png'), 
                   dpi=200, bbox_inches='tight')
        plt.close()

def process_data(base_path, drug_name, user_name):
    """
    Process all concentrations and create feature-based plots
    """
    # Get all concentration folders
    #concentrations = [d for d in os.listdir(base_path) 
    #                 if os.path.isdir(os.path.join(base_path, d)) and d.isdigit()]
    
    # Get available concentrations
    folder_names = [d for d in os.listdir(base_path) if os.path.isdir(os.path.join(base_path, d))]

    # Create a mapping: {original_name: float_value}, ensuring precision is preserved
    concentration_map = {d: float(d) for d in folder_names if _is_valid_float(d)}

    # Sort by numerical value but keep original names
    concentrations = sorted(concentration_map.keys(), key=lambda x: concentration_map[x])

    print("Sorted concentrations:", concentrations)
    print(drug_name)
    print(user_name)
    
    # Read and combine data from all concentrations
    all_data = []
    for concentration in sorted(concentrations, key=float):
        print(f"Processing concentration: {concentration}")
        data = read_feature_data(base_path, drug_name, concentration)
        if data is not None:
            all_data.append(data)
    
    if all_data:
        # Combine all concentration data
        combined_data = pd.concat(all_data, ignore_index=True)
        
        # Create output folder
        output_folder = "./results/plots/features/"
        
        # Create plots
        create_feature_plots(combined_data, drug_name, output_folder)
        print("Created all feature plots")

if __name__ == "__main__":
    # Replace these with your actual values
    if len(sys.argv) < 4: print("Please provide the result folder, drug name and the user name!")
    else:
      print("TEST: " + sys.argv[1])
      temp_string_list = os.path.basename(os.path.dirname(sys.argv[1])).split("_")
      drug_name = sys.argv[2]
      user_name = sys.argv[3]
      print(drug_name)
      print(user_name)
      process_data(sys.argv[1], drug_name, user_name)
