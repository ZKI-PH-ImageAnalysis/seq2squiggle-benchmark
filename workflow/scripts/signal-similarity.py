#!/usr/bin/env python

import logging
import rich_click as click
import sys
import polars as pl
import pyslow5
import numpy as np
from tslearn.metrics import SoftDTWLossPyTorch
from tslearn.metrics import dtw as ts_dtw
from scipy.spatial.distance import euclidean
from fastdtw import fastdtw
import pandas as pd
import pyts.metrics 
from dtaidistance import dtw
from pathlib import Path
import seaborn as sns
import matplotlib.pyplot as plt
import linmdtw
import warnings
warnings.filterwarnings("ignore")
logger = logging.getLogger("signal-similarity")


def path2id(path):
    path_parts = Path(path).with_suffix('').with_suffix('').parts
    id = path_parts[-2] + "-" + path_parts[-1]
    return id



def scatter_dtw(df, ax):
    sns.set_theme(style="whitegrid")
    sns.set_context("paper", font_scale=1.5)
    
    # Determine colors based on position relative to the x=y line
    df['color'] = ['blue' if x < y else 'orange' for x, y in zip(df['DTW_norm_seq2squiggle'], df['DTW_norm_squigulator'])]

    # Create the scatter plot with custom colors
    g = sns.scatterplot(data=df, x="DTW_norm_seq2squiggle", y="DTW_norm_squigulator", hue='color', palette={'blue': '#1f77b4', 'orange': '#ff7f0e'}, legend="full", s=20, ax=ax)
    
    # Set axis labels
    g.set(xlabel='Normalized DTW Deviation (seq2squiggle vs. real)', ylabel='Normalized DTW Deviation (squigulator vs. real)')
    # g.set_title('Comparison of DTW Deviations')
    
    # Draw a line of x=y
    x0, x1 = g.get_xlim()
    y0, y1 = g.get_ylim()
    lims = [max(x0, y0), min(x1, y1)]
    g.plot(lims, lims, '-r')
    
    # Customize legend
    handles, labels = g.get_legend_handles_labels()
    new_labels = ['seq2squiggle better', 'squigulator better']
    g.legend(handles=handles, labels=new_labels, title=None, loc='upper left')

def violin_dtw(df, ax):
    # Prepare data for violin plot
    violin_df = df[["DTW_norm_seq2squiggle", "DTW_norm_squigulator"]].rename(columns={"DTW_norm_seq2squiggle":"seq2squiggle", "DTW_norm_squigulator":"squigulator"})
    violin_df = violin_df.melt(var_name='Tools', value_name='DTW value')
    
    # Plot violin plot
    sns.violinplot(ax=ax, x="Tools", y="DTW value", data=violin_df, hue="Tools")
    ax.set(ylabel='Normalized DTW Deviation')
    #ax.set_title('Comparison of DTW Deviations')

def hist_read_len(df, out):
    # plot bar plot read length
    sns.set_theme(style="whitegrid")
    sns.set_context("paper", font_scale = 1)
    read_len_df = df[["Sim_len_seq2squiggle", "Sim_len_squigulator", "real_len_seq2squiggle"]]
    read_len_df = read_len_df.rename(columns={"Sim_len_seq2squiggle": "Read length seq2squiggle", 
                                "Sim_len_squigulator":"Read length squigulator",
                                "real_len_seq2squiggle":"Read length experimental data"})
    read_len_df = read_len_df.melt(var_name='Tools', value_name='read length')
    hist_df = read_len_df[read_len_df['Tools'] == "Read length experimental data"]
    hist_df2 = read_len_df[read_len_df['Tools'] != "Read length experimental data"]
    ax = sns.displot(read_len_df, x="read length", hue="Tools", multiple="dodge", kde=True, label='small')
    plt.xticks(rotation=45)
    plt.savefig(out, dpi=300)
    plt.close()



@click.group(context_settings=dict(help_option_names=["-h", "--help"]))
def main():
    """
    # signal-similarity

    Calculates the signal similarity between multiple simulated slow5 files and a reference

    SLOW5 needs to be path to a existing slow5 file
    """

@main.command()
@click.option(
    "--seq2squiggle",
    required=True,
    type=click.Path(exists=True, dir_okay=False),
)
@click.option(
    "--squigulator",
    required=True,
    type=click.Path(exists=True, dir_okay=False),
)
@click.option(
    "--out",
    required=True,
    type=click.Path(dir_okay=False),
    help="Path to reference exported picture",
)
def plot(
    seq2squiggle,
    squigulator,
    out
):  
    
    seq2squiggle_df = pd.read_table(seq2squiggle)
    squigulator_df = pd.read_table(squigulator)
    merged_df = seq2squiggle_df.merge(squigulator_df, on="read_id", suffixes=('_seq2squiggle', '_squigulator'))

    median_out = out.replace('.png', '.tsv')
    # Compute median values
    median_values = {
        'median_seq2squiggle': merged_df['DTW_norm_seq2squiggle'].median(),  # Replace 'some_column_seq2squiggle' with your column name
        'median_squigulator': merged_df['DTW_norm_squigulator'].median(),    # Replace 'some_column_squigulator' with your column name
    }
    # Convert median values to DataFrame
    median_df = pd.DataFrame([median_values])
    # Save the median values to TSV
    median_df.to_csv(median_out, sep='\t', index=False)



    sns.set_theme(style="whitegrid")
    sns.set_context("paper", font_scale=1.5)

    # Create figure and axes
    fig, axes = plt.subplots(1, 2, figsize=(15, 6))

    # Plot scatter plot on the left axis
    scatter_dtw(merged_df, axes[0])
    
    # Plot violin plot on the right axis
    violin_dtw(merged_df, axes[1])
    
    # Add subplot labels (A and B)
    axes[0].text(-0.1, 1.04, 'A', transform=axes[0].transAxes, size=20, weight='bold')
    axes[1].text(-0.1, 1.04, 'B', transform=axes[1].transAxes, size=20, weight='bold')
    
    # Remove subplot titles
    axes[0].set_title('')
    axes[1].set_title('')
    
    # Adjust layout and save the figure
    plt.tight_layout()
    plt.savefig(out, dpi=600)
    plt.close()



@main.command()
@click.argument(
    "slow5",
    required=True,
    nargs=1,
    type=click.Path(exists=True, dir_okay=False),
)
@click.option(
    "--reference",
    required=True,
    type=click.Path(exists=True, dir_okay=False),
    help="Path to reference signal file",
)
@click.option(
    "--out",
    required=True,
    type=click.Path(dir_okay=False),
    help="Path to reference output tsv",
)
@click.option(
    "--dtw_impl",
    type=click.Choice(['fastdtw', 'dtw-tsl', 'dtw-dtaid', 'dtw-pyts', 'linmdtw-fastdtw', 'linmdtw-mrmsdtw']), 
    default="linmdtw-fastdtw",
    help="Type of dtw implementation to use.",
)
@click.option(
    "--head_n",
    type=int,
    default=5000,
    help="How many samples to process from the joined df"
)
@click.option(
    "--norm",
    type=bool,
    default=True,
    help="If dtw value should be normalized by dividing through max length of signal"
)
def compare(
    slow5,
    reference,
    out,
    dtw_impl,
    head_n,
    norm
):
    
    setup_logging("info")
    read_list = []
    dtw_list = []
    # Open slow5 simulated file
    s5_sim = pyslow5.Open(slow5,'r')
    s5_ref = pyslow5.Open(reference, 'r')
    reads_sim = s5_sim.seq_reads()
    for idx, read_sim in enumerate(reads_sim):
        signal_sim = read_sim['range']/read_sim['digitisation'] * (read_sim['signal'] + read_sim['offset'])
        read_ref = s5_ref.get_read(read_sim['read_id'], aux='all')
        if read_ref == None:
            continue
        signal_ref = read_ref['range']/read_ref['digitisation'] * (read_ref['signal'] + read_ref['offset'])

        if dtw_impl == "fastdtw":
            signal_sim = np.asarray(signal_sim).flatten().reshape((-1, 1))
            signal_ref = np.asarray(signal_ref).flatten().reshape((-1, 1))     
            dtw_dist, _ = fastdtw(signal_sim, signal_ref, dist=euclidean)
        elif dtw_impl == "dtw-pyts":
            signal_sim = np.asarray(signal_sim).flatten()
            signal_ref = np.asarray(signal_ref).flatten()
            dtw_dist = pyts.metrics.dtw(signal_sim, signal_ref, method='fast', return_path=False)
        elif dtw_impl == "dtw-dtaid":
            signal_sim = np.asarray(signal_sim, dtype=np.double).flatten()
            signal_ref = np.asarray(signal_ref, dtype=np.double).flatten()
            dtw_dist = dtw.distance_fast(signal_sim, signal_ref)
        elif dtw_impl == "dtw-tsl":
            signal_sim = np.asarray(signal_sim).flatten().reshape((-1, 1))
            signal_ref = np.asarray(signal_ref).flatten().reshape((-1, 1))
            dtw_dist = ts_dtw(signal_sim, signal_ref)
        elif dtw_impl == "linmdtw-fastdtw":
            signal_sim = np.asarray(signal_sim).flatten().reshape((-1, 1))
            signal_ref = np.asarray(signal_ref).flatten().reshape((-1, 1))
            dtw_dist, _ = linmdtw.fastdtw(signal_sim, signal_ref, radius=5)
        elif dtw_impl == "linmdtw-mrmsdtw":
            signal_sim = np.asarray(signal_sim).flatten().reshape((-1, 1))
            signal_ref = np.asarray(signal_ref).flatten().reshape((-1, 1))
            dtw_dist, _ = linmdtw.mrmsdtw(signal_sim, signal_ref, tau=10**7)
            # dtw_dist, _ = linmdtw.mrmsdtw(signal_sim, signal_ref, tau=10**4)

        if norm:
            max_len = max(len(signal_sim), len(signal_ref))
            dtw_dist = dtw_dist / max_len  

        read_list.append(read_sim['read_id'])
        dtw_list.append(dtw_dist)

        if idx == head_n:
            break
    s5_sim.close()
    s5_ref.close()

    data = {'read_id': read_list,
            'DTW_norm': dtw_list,
    }
    df_out = pl.DataFrame._from_dict(data)
    df_out.write_csv(out, separator="\t")

def setup_logging(verbosity):
    logging_levels = {
        "debug": logging.DEBUG,
        "info": logging.INFO,
        "warning": logging.WARNING,
        "error": logging.ERROR,
    }

    # Configure logging.
    logging.captureWarnings(True)
    root_logger = logging.getLogger()
    root_logger.setLevel(logging.DEBUG)
    warnings_logger = logging.getLogger("py.warnings")

    # Formatters for file vs console:
    console_formatter = logging.Formatter(
        "{name} {levelname} {asctime}: {message}", style="{", datefmt="%H:%M:%S"
    )

    console_handler = logging.StreamHandler(sys.stderr)
    console_handler.setLevel(logging_levels[verbosity.lower()])
    console_handler.setFormatter(console_formatter)
    root_logger.addHandler(console_handler)
    warnings_logger.addHandler(console_handler)

    logging.getLogger("fsspec").setLevel(logging.WARNING)
    logging.getLogger("github").setLevel(logging.WARNING)
    logging.getLogger("h5py").setLevel(logging.WARNING)
    logging.getLogger("numba").setLevel(logging.WARNING)
    logging.getLogger("pytorch_lightning").setLevel(logging.WARNING)
    logging.getLogger("torch").setLevel(logging.WARNING)
    logging.getLogger("urllib3").setLevel(logging.WARNING)

if __name__ == "__main__":
    main()