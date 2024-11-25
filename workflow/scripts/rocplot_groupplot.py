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

def plot_grouped_figure(data_r9, tools_r9, data_r10, tools_r10, output_path):
    sns.set_theme(style="whitegrid")
    sns.set_context("paper", font_scale=1.5)
    
    # Define custom color palette to skip the third color
    custom_palette = sns.color_palette()
    selected_colors_r9 = [custom_palette[0], custom_palette[1], custom_palette[3], custom_palette[4]]
    selected_colors_r10 = [custom_palette[0], custom_palette[1]]

    fig, axes = plt.subplots(2, 2, figsize=(14, 12))

    # R9 plots
    for idx, (data, tool_name) in enumerate(zip(data_r9, tools_r9)):
        data = data[(data['sensitivity'] != 0.0) | (data['precision'] != 0.0)]
        
        # True Positives vs False Positives (Top Left - A)
        sns.lineplot(ax=axes[0, 0], x=data['false_positives'], y=data['true_positives_baseline'],
                     label=tool_name, lw=2.0, color=selected_colors_r9[idx])
        
        # Precision-Recall Curve (Top Right - B)
        sns.lineplot(ax=axes[0, 1], x=data['sensitivity'], y=data['precision'],
                     label=tool_name, lw=2.0, color=selected_colors_r9[idx])

    # R10 plots
    for idx, (data, tool_name) in enumerate(zip(data_r10, tools_r10)):
        data = data[(data['sensitivity'] != 0.0) | (data['precision'] != 0.0)]
        
        # True Positives vs False Positives (Bottom Left - C)
        sns.lineplot(ax=axes[1, 0], x=data['false_positives'], y=data['true_positives_baseline'],
                     label=tool_name, lw=2.0, color=selected_colors_r10[idx])
        
        # Precision-Recall Curve (Bottom Right - D)
        sns.lineplot(ax=axes[1, 1], x=data['sensitivity'], y=data['precision'],
                     label=tool_name, lw=2.0, color=selected_colors_r10[idx])

    # Customize each axis
    axes[0, 0].set_xlim(left=0.0)
    axes[0, 0].set_ylim(bottom=0.0)
    axes[0, 0].set_xlabel('False Positives')
    axes[0, 0].set_ylabel('True Positives')
    axes[0, 0].set_title('True Positives vs. False Positives (R9.4.1)')
    axes[0, 0].xaxis.set_major_formatter(FuncFormatter(format_large_numbers))
    axes[0, 0].yaxis.set_major_formatter(FuncFormatter(format_large_numbers))
    axes[0, 0].text(-0.1, 1.1, 'A', transform=axes[0, 0].transAxes, size=20, weight='bold')

    axes[0, 1].set_xlim([0.0, 1.05])
    axes[0, 1].set_ylim([0.0, 1.05])
    axes[0, 1].set_ylabel('Precision')
    axes[0, 1].set_xlabel('Recall')
    axes[0, 1].set_title('Precision-Recall Curve (R9.4.1)')
    axes[0, 1].text(-0.1, 1.1, 'B', transform=axes[0, 1].transAxes, size=20, weight='bold')

    # Set R9.4.1 title for A + B
    #axes[0, 0].set_title("R9.4.1", loc='left', size=16, weight='bold', pad=20)
    #axes[0, 1].set_title("R9.4.1", loc='left', size=16, weight='bold', pad=20)

    axes[1, 0].set_xlim(left=0.0)
    axes[1, 0].set_ylim(bottom=0.0)
    axes[1, 0].set_xlabel('False Positives')
    axes[1, 0].set_ylabel('True Positives')
    axes[1, 0].set_title('True Positives vs. False Positives (R10.4.1)')
    axes[1, 0].xaxis.set_major_formatter(FuncFormatter(format_large_numbers))
    axes[1, 0].yaxis.set_major_formatter(FuncFormatter(format_large_numbers))
    axes[1, 0].text(-0.1, 1.1, 'C', transform=axes[1, 0].transAxes, size=20, weight='bold')

    axes[1, 1].set_xlim([0.0, 1.05])
    axes[1, 1].set_ylim([0.0, 1.05])
    axes[1, 1].set_ylabel('Precision')
    axes[1, 1].set_xlabel('Recall')
    axes[1, 1].set_title('Precision-Recall Curve (R10.4.1)')
    axes[1, 1].text(-0.1, 1.1, 'D', transform=axes[1, 1].transAxes, size=20, weight='bold')

    # Set R10.4.1 title for A + B
    #axes[1, 0].set_title("R10.4.1", loc='left', size=16, weight='bold', pad=20)
    #axes[1, 1].set_title("R10.4.1", loc='left', size=16, weight='bold', pad=20)

    # Create the shared legend with the title "Tools"
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc='lower center', ncol=4, frameon=True, bbox_to_anchor=(0.5, 0.01))

    # Remove legends from individual subplots
    for ax in axes.flatten():
        ax.get_legend().remove()

    plt.tight_layout(rect=[0, 0.07, 1, 1])  
    plt.subplots_adjust(hspace=0.4)

    plt.savefig(output_path, dpi=600)
    plt.close()

def main():
    parser = argparse.ArgumentParser(description="Generate grouped figure for R9 and R10 data.")
    # R9 data and tools
    parser.add_argument('--r9_files', nargs='+', required=True, 
                        help="Paths to the TSV files for R9 data (space-separated).")
    parser.add_argument('--r9_tools', nargs='+', required=True,
                        help="Names of tools corresponding to each R9 file (space-separated).")
    # R10 data and tools
    parser.add_argument('--r10_files', nargs='+', required=True,
                        help="Paths to the TSV files for R10 data (space-separated).")
    parser.add_argument('--r10_tools', nargs='+', required=True,
                        help="Names of tools corresponding to each R10 file (space-separated).")
    # Output path
    parser.add_argument('--outdir', required=True, help="Output path for the saved plot.")

    args = parser.parse_args()

    # Ensure output directory exists
    os.makedirs(os.path.dirname(args.outdir), exist_ok=True)

    # Load R9 data using load_data function
    data_r9 = []
    for file_path in args.r9_files:
        data = load_data(file_path)
        if data is not None:
            data_r9.append(data)
        else:
            print(f"Skipping R9 file '{file_path}' due to loading error.")
            return

    # Load R10 data using load_data function
    data_r10 = []
    for file_path in args.r10_files:
        data = load_data(file_path)
        if data is not None:
            data_r10.append(data)
        else:
            print(f"Skipping R10 file '{file_path}' due to loading error.")
            return

    # Validate file and tool count matches
    if len(data_r9) != len(args.r9_tools):
        print("Error: The number of R9 files does not match the number of R9 tool names provided.")
        return
    if len(data_r10) != len(args.r10_tools):
        print("Error: The number of R10 files does not match the number of R10 tool names provided.")
        return

    # Plot the grouped figure
    plot_grouped_figure(data_r9, args.r9_tools, data_r10, args.r10_tools, args.outdir)

if __name__ == "__main__":
    main()