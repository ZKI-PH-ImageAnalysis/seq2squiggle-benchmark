import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import argparse

def read_data(input_file):
    # Read data from file
    data = pd.read_csv(input_file, sep='\t')
    return data

def bytes_to_megabytes(bytes_value):
    return bytes_value / (1024 * 1024)

def seconds_to_minutes(seconds_value):
    return seconds_value / 60.0

def plot_comparison(squigulator, seq2squiggle, outfile):
    # Read data from both files
    data1 = read_data(squigulator)
    data2 = read_data(seq2squiggle)
    file_names = ["squigulator", "seq2squiggle"]

    # Melt dataframes to long format
    data1_melt = data1.melt(id_vars=['h:m:s'], var_name='Metric', value_name='Value')
    data2_melt = data2.melt(id_vars=['h:m:s'], var_name='Metric', value_name='Value')

    # Convert runtime values to minutes
    data1_melt.loc[data1_melt['Metric'] == 's', 'Value'] = data1_melt[data1_melt['Metric'] == 's']['Value'].apply(seconds_to_minutes)
    data2_melt.loc[data2_melt['Metric'] == 's', 'Value'] = data2_melt[data2_melt['Metric'] == 's']['Value'].apply(seconds_to_minutes)

    # Convert memory values to megabytes
    data1_melt.loc[data1_melt['Metric'] == 'max_rss', 'Value'] = data1_melt[data1_melt['Metric'] == 'max_rss']['Value'].apply(bytes_to_megabytes)
    data2_melt.loc[data2_melt['Metric'] == 'max_rss', 'Value'] = data2_melt[data2_melt['Metric'] == 'max_rss']['Value'].apply(bytes_to_megabytes)

    # Combine dataframes
    combined_data = pd.concat([data1_melt, data2_melt], keys=file_names)

    # Plot using Seaborn
    plt.figure(figsize=(12, 8))

    # Runtime comparison
    plt.subplot(2, 1, 1)
    sns.barplot(data=combined_data[combined_data['Metric'] == 's'], x=file_names, y='Value', hue=file_names)
    plt.xlabel('Iterations')
    plt.ylabel('Runtime (minutes)')
    plt.title('Runtime Comparison')

    # Memory consumption comparison
    plt.subplot(2, 1, 2)
    sns.barplot(data=combined_data[combined_data['Metric'] == 'max_rss'], x=file_names, y='Value', hue=file_names)
    plt.xlabel('Tools')
    plt.ylabel('Max RSS (in Megabytes)')
    plt.title('Max RSS Comparison')

    plt.tight_layout()
    plt.savefig(outfile, dpi=600)



if __name__ == "__main__":
    # Parse command-line arguments
    parser = argparse.ArgumentParser(description='Plot runtime and memory consumption comparison.')
    parser.add_argument('squigulator', type=str, help='Path to squigulator file')
    parser.add_argument('seq2squiggle', type=str, help='Path to seq2squiggle file')
    parser.add_argument('outfile', type=str, help='Path to png file')
    args = parser.parse_args()

    # Plot comparison
    plot_comparison(args.squigulator, args.seq2squiggle, args.outfile)
