#!/usr/bin/env python3

import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import random
import sys
from pathlib import Path

# Define a global variable
random_sample = []
fixed_sample = [1,2,3,4,5]

def generate_random_sample(number_of_data):
    global random_sample  # Declare the global variable
    random_sample = random.sample(range(0, number_of_data), 5)


# Helper function to check if a string can be converted to float
def _is_valid_float(s):
    try:
        float(s)
        return True
    except ValueError:
        return False

def read_concentration_data(folder_path, drug_name, concentration, num_samples=10, is_control=False):
    """
    Read data for a specific concentration from multiple sample files.
    """
    data_dict = {}
    print(concentration)

    for i in random_sample:
        if is_control:
            i = 0  # Control sample is always 0
        pattern = "{}_{:.2f}_time_series_smp{}.csv".format(drug_name, float(concentration), i)
        print(pattern)
        print(folder_path)
        matching_files = list(Path(folder_path).glob(pattern))
        print(matching_files)
        if matching_files:
            file_path = matching_files[0]  # Take the first matching file
            try:
                df = pd.read_csv(file_path)
                data_dict[f'smp{i}'] = df

            except Exception as e:
                print(f"Error reading file {file_path}: {e}")
        
        if is_control:
            break  # Only need one control sample

    return data_dict

def create_multi_panel_plot(base_path, drug_name, concentrations, feature_name, ylabel, title, sample_size):
    """
    Create a multi-panel plot for a specific feature across different concentrations.
    """
    # Sort concentrations to ensure correct panel placement
    #concentrations = sorted([int(c) for c in concentrations if int(c) > 0])
    concentrations = sorted([float(c) for c in concentrations if float(c) > 0])

    if not concentrations:
        print("No valid concentrations found")
        return

    # Calculate rows and columns for the subplots
    n_panels = len(concentrations)
    n_rows = (n_panels + 1) // 2
    n_cols = 2
    
    # Create figure
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(8, 6))
    #fig.suptitle(f'{title} across different concentrations')
    
    # Flatten axes array for easier iteration
    axes_flat = axes.flatten()
    
    # Find any folder that represents "0" in float form
    control_folder = conc_str = f"0.00"  # Ensure two decimal places

    if control_folder:
      control = read_concentration_data(os.path.join(base_path, control_folder), drug_name, 0, num_samples=1, is_control=True)
    else:
      print("No control folder found.")
      control = {}

    
    # Plot each concentration in its panel
    for idx, conc in enumerate(concentrations):
        ax = axes_flat[idx]
        conc_str = f"{conc:.2f}"  # Ensure two decimal places
        data_dict = read_concentration_data(os.path.join(base_path, conc_str), drug_name, conc, sample_size)
        
        if data_dict:
            # Plot control first (black dashed line)
            for c_sample, c_df in control.items():
                ax.plot(c_df['Time(ms)'], c_df[feature_name], 
                       'k--', label='Control', alpha=0.7, linewidth=1)
            
            # Plot samples
            for sample, df in data_dict.items():
                ax.plot(df['Time(ms)'], df[feature_name], 
                       alpha=0.7, linewidth=1, label=f'Sample {sample.split("smp")[1]}')
            
            ax.set_title(f'Concentration: {conc}')
            ax.set_xlabel('Time (ms)')
            ax.set_ylabel(ylabel)
            ax.grid(True, alpha=0.3)
    
    # Remove empty subplots if any
    for idx in range(len(concentrations), len(axes_flat)):
        fig.delaxes(axes_flat[idx])
    
    # Add legend to the figure
    # Get lines and labels from first subplot
    lines = axes_flat[0].get_lines()
    labels = ['Control'] + [f'Sample {i}' for i in fixed_sample]
    
    # Place legend in the center of the figure
    fig.legend(lines[:6], labels, 
              loc='center', 
              bbox_to_anchor=(0.5, 0.5))
    
    plt.tight_layout()
    return fig

def plot_all_features(base_path, drug_name, sample_size):
    """
    Plot all features from the data files.
    """
    # Get available concentrations
    #concentrations = [d for d in os.listdir(base_path) 
    #                 if os.path.isdir(os.path.join(base_path, d)) and d.isdigit()]
    #print("concs: ", concentrations)

    # Get available concentrations
    folder_names = [d for d in os.listdir(base_path) if os.path.isdir(os.path.join(base_path, d))]

    # Create a mapping: {original_name: float_value}, ensuring precision is preserved
    concentration_map = {d: float(d) for d in folder_names if _is_valid_float(d)}

    # Sort by numerical value but keep original names
    sorted_concentrations = sorted(concentration_map.keys(), key=lambda x: concentration_map[x])

    print("Sorted concentrations:", sorted_concentrations)

    
    # Dictionary of features and their labels
    features = {
        'Vm(mV)': 'Voltage (mV)',
        'dVm/dt(mV/ms)': 'dV/dt (mV/ms)',
        'Cai(nM)': 'Cai (nM)',
        'INa(nA)': 'INa (nA)',
        'INaL(nA)': 'INaL (nA)',
        'ICaL(nA)': 'ICaL (nA)',
        'Ito(nA)': 'Ito (nA)',
        'IKr(nA)': 'IKr (nA)',
        'IKs(nA)': 'IKs (nA)',
        'IK1(nA)': 'IK1 (nA)',
        'Inet(uA)': 'iNet (uA)',
        'Inet_AUC(uC)': 'AUC of iNet (uC)'
    }
    
    # Create output directory if it doesn't exist
    print("TEST0  " + base_path)
    print("TEST1  " + drug_name)
    output_folder_name = "./results/plots/time_series/"
    print("TEST2  " + output_folder_name)
    os.makedirs(output_folder_name, exist_ok=True)
    
    # Create a plot for each feature
    for feature, ylabel in features.items():
        fig = create_multi_panel_plot(base_path, drug_name, sorted_concentrations, feature, 
                                    ylabel, feature.split('(')[0], sample_size)
        if fig:
            file_name = f"{feature.replace('/','_').split('(')[0]}_plot.png"
            plt_fullname = output_folder_name + file_name
            plt.savefig(plt_fullname, dpi=200, bbox_inches='tight')
            plt.close(fig)

# Example usage
if __name__ == "__main__":
    # Replace with your actual base path
    if len(sys.argv) < 4: print("Please provide the result folder path, drug name and number of sample!!!")
    else:
      if (int(sys.argv[3]) <= 0): print("Incorrect data size --- {}!!".format(sys.argv[3]))
      else:
        generate_random_sample(int(sys.argv[3]))
        plot_all_features(sys.argv[1], sys.argv[2], sys.argv[3])
