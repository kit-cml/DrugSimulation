import os
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

def scatter_plots_from_df(df,plt_file,input_folder):
    # Exclude non-numeric columns (if any) and 'Sample_id'
    numeric_cols = df.select_dtypes(include='number').columns
    
    # Create scatter plots for each numeric column
    for col in numeric_cols:
        plt.figure()
        plt.scatter(df['Sample_ID'], df[col])
        plt.title("Scatter plot of {} {} {}".format(col, plt_file.split("_")[0], plt_file.split("_")[1]))
        plt.xlabel('Sample ID')
        plt.ylabel(col)
        plt.gca().xaxis.set_major_formatter(ticker.FuncFormatter(lambda x, _: '{:.0f}'.format(x)))
        plt.grid(False)
        plt.savefig("./plots/"+input_folder.split("/")[2]+"/"+input_folder.split("/")[3]+"/"+col+".png")
        plt.close()


def read_and_sort_plt_files(input_folder):
    # Get a list of all .plt files in the input folder
    plt_files = [f for f in os.listdir(input_folder) if f.endswith('.plt') and "features" in f]
    
    # Initialize an empty list to store DataFrames
    df_list = []
    
    # Process each .plt file
    for plt_file in plt_files:
        # Read data from .plt file using pandas
        df = pd.read_csv(os.path.join(input_folder, plt_file))
        df_list.append(df)
    
    # Concatenate all DataFrames into a single DataFrame
    combined_df = pd.concat(df_list, ignore_index=True)
    
    # Sort DataFrame based on "Sample_id" column
    combined_df_sorted = combined_df.sort_values(by='Sample_ID')
    scatter_plots_from_df(combined_df_sorted, plt_file, input_folder)
    combined_df_sorted.to_csv(input_folder+"/feature_joined_"+plt_file.split("_")[0]+"_"+plt_file.split("_")[1]+".plt", index=False)

# Example usage:
input_folder = './result/'
for drugname in os.listdir(input_folder):
	for conc in os.listdir(os.path.join(input_folder,drugname)):
		read_and_sort_plt_files(os.path.join(input_folder,drugname,conc))
