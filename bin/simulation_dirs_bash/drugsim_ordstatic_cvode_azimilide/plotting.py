import os
import matplotlib.pyplot as plt
import pandas as pd

def plot_from_plt_files(input_folder, output_folder, fig, axs, col_num):
    # Ensure the output folder exists
	os.makedirs(output_folder, exist_ok=True)
    
    # Get a list of all .plt files in the input folder
	plt_files = [f for f in os.listdir(input_folder) if f.endswith('.plt') and "time_series" in f]
    # Process each .plt file
	for plt_file in plt_files:
        # Read data from .plt file (assuming it's a text file with x, y values)
		print(plt_file)
		df = pd.read_csv(os.path.join(input_folder, plt_file))
		for idx, col in enumerate(df.columns, start=0):
			if col != 'Time':
				print(col)
				idx = idx-1
				print(idx)
	        		# Create plot
				#plt.figure()
				axs[idx,col_num].plot(df["Time"], df[col])
				#axs[idx,col_num].title('{} of {} at {} mMol'.format(col.split("(")[0],plt_file.split("_")[0], plt_file.split("_")[1]))
				axs[idx,col_num].set(xlabel ='time (msec)', ylabel = col)
				#axs[idx,col_num].ylabel(col)
				#plt.legend()
				if "/" in col: col=col.replace("/","")
				#plot_name = os.path.splitext(plt_file)[0] + '_plot_' + col.split(" ")[0] + '.png'  # Change extension to .png
				#plt.savefig(os.path.join(output_folder, plot_name))
				#plt.close()


def all_plot_from_plt_files(input_folder, output_folder):
    # Ensure the output folder exists
        os.makedirs(output_folder, exist_ok=True)

    # Get a list of all .plt files in the input folder
        plt_files = [f for f in os.listdir(input_folder) if f.endswith('.plt') and "time_series" in f]

    # Process each .plt file
        for plt_file in plt_files:
        # Read data from .plt file (assuming it's a text file with x, y values)
                #print(plt_file)
                df = pd.read_csv(os.path.join(input_folder, plt_file))
                for col in df.columns:
                        if col != 'Time':
                               # print(col)
                                # Create plot
                                # plt.figure()
                                plt.plot(df["Time"], df[col])
                                #plt.title('Plot from {} of {} @ {} mMol'.format(col.split("(")[0],plt_file.split("_")[0], plt_file.split("_")[1]))
                                plt.xlabel('time (msec)')
                                plt.ylabel(col)
                                plt.legend()   	
                                if "/" in col: col=col.replace("/","")
                                #plot_name = os.path.splitext(plt_file)[0] + '_plot_' + col.split(" ")[0] + '.png'  # Change extension to .png
                                #plt.savefig(os.path.join(output_folder, plot_name))
                                #plt.close()



# Example usage:
input_root = './result'
output_root = './plots'
col_num = 0
for drugname in os.listdir(input_root):
	fig, axs = plt.subplots(14, 10, figsize=(16, 30))
	for conc in os.listdir(os.path.join(input_root, drugname)):
		print(conc)
		joined = os.path.join(drugname,conc)
		input_folder = os.path.join(input_root, joined)
		output_folder = os.path.join(output_root, joined)
		plot_from_plt_files(input_folder, output_folder, fig, axs, col_num)
		col_num = col_num + 1
	plot_name = drugname + '.png'  # Change extension to .png
	col_num = 0
	fig.savefig(os.path.join(output_root, plot_name))
	print(os.path.join(output_root,plot_name))
	#fig.close()
