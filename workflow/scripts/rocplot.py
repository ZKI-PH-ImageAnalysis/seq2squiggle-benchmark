import argparse
import gzip
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.metrics import roc_curve, auc
import os
from matplotlib.ticker import FuncFormatter

def load_data(file_path):
    with gzip.open(file_path, 'rt') as f:
        # Read all lines
        lines = f.readlines()
        
        # Find the last comment line starting with '#'
        header_line = None
        for line in lines:
            if line.startswith('#') and not line.startswith('#Version'):
                header_line = line
        if header_line is None:
            raise ValueError("Header line not found")

        # Extract the header from the last comment line
        header = header_line.lstrip('#').strip().split()
        
        # Read the data into a DataFrame
        data = pd.read_csv(file_path, sep='\t', comment='#', names=header, skiprows=1)
    
    return data

def plot_roc(data1, data2, output_path):
    # Compute TPR and FPR
    tpr1 = data1['true_positives_baseline'] / (data1['true_positives_baseline'] + data1['false_negatives'])
    fpr1 = data1['false_positives'] / (data1['false_positives'].max() + data1['true_positives_baseline'].max())
    roc_auc1 = auc(fpr1, tpr1)

    tpr2 = data2['true_positives_baseline'] / (data2['true_positives_baseline'] + data2['false_negatives'])
    fpr2 = data2['false_positives'] / (data2['false_positives'].max() + data2['true_positives_baseline'].max())
    roc_auc2 = auc(fpr2, tpr2)

    # Plotting with Seaborn
    sns.set_theme(style="whitegrid")
    sns.set_context("paper", font_scale=1.0)

    plt.figure(figsize=(6, 6))
    sns.lineplot(x=fpr1, y=tpr1, label='seq2squiggle (area = %0.2f)' % roc_auc1, lw=1.5)
    sns.lineplot(x=fpr2, y=tpr2, label='squigulator (area = %0.2f)' % roc_auc2, lw=1.5)
    #plt.plot([0, 1], [0, 1], color='gray', lw=2, linestyle='--')

    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title('Receiver Operating Characteristic')
    plt.legend(loc="lower right")
    plt.tight_layout()
    plt.savefig(output_path, dpi=600)
    plt.close()

def plot_specificity_sensitivity(data1, data2, output_path):
    # Compute Sensitivity and Specificity
    sns.set_theme(style="whitegrid")
    sns.set_context("paper", font_scale=1.0)

    plt.figure(figsize=(6, 6))
    sns.lineplot(x=data1['false_positives'], y=data1['true_positives_baseline'], label='seq2squiggle', lw=1.5)
    sns.lineplot(x=data2['false_positives'], y=data2['true_positives_baseline'], label='squigulator', lw=1.5)

    plt.xlim([0.0, 1.05])
    plt.ylim([0.0, 1.05])
    plt.ylabel('Specificity')
    plt.xlabel('Sensitivity')
    plt.title('Specificity-Sensitivity Curve')
    plt.legend(loc="lower left")
    plt.tight_layout()
    plt.savefig(output_path, dpi=600)
    plt.close()

def plot_true_false_positives(data1, data2, output_path):
    sns.set_theme(style="whitegrid")
    sns.set_context("paper", font_scale=1.0)

    plt.figure(figsize=(6, 6))
    sns.lineplot(x=data1['false_positives'], y=data1['true_positives_baseline'], label='seq2squiggle', lw=1.5)
    sns.lineplot(x=data2['false_positives'], y=data2['true_positives_baseline'], label='squigulator', lw=1.5)

    plt.xlim(left=0.0)
    plt.ylim(bottom=0.0)

    plt.xlabel('False Positives')
    plt.ylabel('True Positives')
    plt.title('True Positives vs. False Positives')
    plt.legend(loc="upper left")
    plt.tight_layout()
    plt.savefig(output_path, dpi=600)
    plt.close()

def format_large_numbers(value, tick_number):
    if value >= 1_000_000:
        return f'{value / 1_000_000:.1f}M'
    elif value >= 1_000:
        return f'{value / 1_000:.1f}K'
    else:
        return str(int(value))

def plot_grouped(data_list, tool_names, output_path):
    sns.set_theme(style="whitegrid")
    sns.set_context("paper", font_scale=1.25)

    # Define color palette, skipping the third default color
    custom_palette = sns.color_palette() # Use colorblind for default colors
    selected_colors = [custom_palette[0], custom_palette[1], custom_palette[3], custom_palette[4]]


    fig, axes = plt.subplots(1, 2, figsize=(12, 6))

    for idx, (data, tool_name) in enumerate(zip(data_list, tool_names)):
        data = data[(data['sensitivity'] != 0.0) | (data['precision'] != 0.0)]

        # True Positives vs. False Positives Plot
        sns.lineplot(ax=axes[0], x=data['false_positives'], y=data['true_positives_baseline'], label=tool_name, lw=2.0, color=selected_colors[idx])

        # Precision-Recall Plot
        sns.lineplot(ax=axes[1], x=data['sensitivity'], y=data['precision'], label=tool_name, lw=2.0, color=selected_colors[idx])
    
    axes[0].set_xlim(left=0.0)
    axes[0].set_ylim(bottom=0.0)
    axes[0].set_xlabel('False Positives')
    axes[0].set_ylabel('True Positives')
    axes[0].set_title('True Positives vs. False Positives')
    axes[0].xaxis.set_major_formatter(FuncFormatter(format_large_numbers))
    axes[0].yaxis.set_major_formatter(FuncFormatter(format_large_numbers))
    axes[0].text(-0.1, 1.1, 'A', transform=axes[0].transAxes, size=20, weight='bold')

    axes[1].set_xlim([0.0, 1.05])
    axes[1].set_ylim([0.0, 1.05])
    axes[1].set_ylabel('Precision')
    axes[1].set_xlabel('Recall')
    axes[1].set_title('Precision-Recall Curve')
    axes[1].text(-0.1, 1.1, 'B', transform=axes[1].transAxes, size=20, weight='bold')

    # Remove individual legends and add a shared legend
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc='lower center', ncol=2, frameon=True, bbox_to_anchor=(0.5, -0.015))
    for ax in axes.flatten():
        ax.get_legend().remove()

    plt.tight_layout(rect=[0, 0.05, 1, 1])
    plt.subplots_adjust(hspace=0.3)
    
    plt.savefig(output_path, dpi=600)
    plt.close()

def main():
    parser = argparse.ArgumentParser(description='Generate ROC, Specificity-Sensitivity, and True/False Positives plots from TSV files.')
    parser.add_argument('--outdir', required=True, help='Output directory for the PNG files')
    parser.add_argument('inputs', nargs='+', help='Input TSV files (gzipped)')
    parser.add_argument('--tools', nargs='+', help='Names of the tools for labeling the datasets')
    args = parser.parse_args()

    # Ensure the output directory exists
    os.makedirs(args.outdir, exist_ok=True)

    # Define the output file paths
    roc_output_path = os.path.join(args.outdir, 'roc_output.png')
    ss_output_path = os.path.join(args.outdir, 'ss_output.png')
    tf_output_path = os.path.join(args.outdir, 'tf_output.png')
    grouped_output_path = os.path.join(args.outdir, 'Figure04.png')

    # Load data for each input file
    data_list = [load_data(input_file) for input_file in args.inputs]

    plot_grouped(data_list, args.tools, grouped_output_path)

if __name__ == "__main__":
    main()