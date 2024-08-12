import polars as pl
from collections import Counter
from multiprocessing import Pool
import matplotlib.pyplot as plt
import argparse
from itertools import product
import seaborn as sns
import os
from matplotlib.collections import PathCollection

def extract_9mers(tsv_file, max_rows=100000000):
    # Read csv
    df = pl.read_csv(tsv_file, separator='\t', n_rows=max_rows)
    # Filter out all NNNNNNNN
    df = df.sort(["position"]).filter(pl.col("model_kmer") != "NNNNNNNNN")
    # Extract the "sequence" column
    sequence_column = df['model_kmer'].to_list()

    return df, sequence_column

def plot_frequency_distribution(kmers, output_file, max_len=10000000):
    kmers = kmers[:max_len]
    # Count occurrences of each 9-mer
    nine_mer_counts = Counter(kmers)
    # Extract unique 9-mers and their counts and sort them
    sorted_counts = sorted(nine_mer_counts.items(), key=lambda x: x[1], reverse=True)
    unique_nine_mers, counts = zip(*sorted_counts)
    max_c = max(counts)
    indices = range(len(unique_nine_mers))

    # Plot the occurrences
    fig = plt.figure(figsize=(10, 6))
    ax = fig.add_subplot(111)
    # sns.barplot(x=indices, y=counts)
    plt.bar(indices, counts)
    plt.plot(indices, counts, color='blue')
    ax.set_yscale('log')
    plt.margins(x=0) 
    plt.xlabel('9-mer')
    plt.ylabel('Occurrences')
    plt.title('Occurrences of 9-mers')
    plt.ylim(0, max_c)
    plt.grid(axis='y')
    plt.xticks(rotation=45)
    plt.tight_layout() 
    plt.savefig(output_file, dpi=600)  # Save the plot to a PNG file
    plt.close()


def plot_boxenplot(
    dataset,
    selected_kmers=["AAAAA", "AATAA", "AACAA", "AAGAA"],
    output_name="subset-violin01",
    y="model_stdv",
):
    kmer_subset01 = dataset[dataset.model_kmer.isin(selected_kmers)]
    fig, ax = plt.subplots(figsize=(24, 8))
    sns.set_theme(style="ticks")
    ax = sns.boxenplot(
        data=kmer_subset01,
        x="model_kmer",
        y=y,
        width=0.25,
        hue="model_kmer",
        palette="turbo",
    )

    plt.minorticks_on()
    plt.ylim(0, 45)
    plt.gca().yaxis.grid(True, which="major")
    plt.gca().yaxis.grid(True, which="minor", linewidth=0.2)
    fig.savefig(output_name, dpi=600)
    


def plot_violin(
    dataset,
    selected_kmers=["AAAAA", "AATAA", "AACAA", "AAGAA"],
    output_name="subset-violin01",
    y="model_stdv",
):
    kmer_subset01 = dataset[dataset.model_kmer.isin(selected_kmers)]
    fig, ax = plt.subplots(figsize=(24, 8))
    sns.set_theme(style="ticks")
    ax = sns.violinplot(
        data=kmer_subset01,
        x="model_kmer",
        y=y,
        width=0.45,
        hue="model_kmer",
        saturation=0.4,
        linewidth=0,
        palette="turbo",
        inner=None,
        cut=0,
        bw_adjust=0.5
    )
    ax = sns.boxplot(
        data=kmer_subset01,
        x="model_kmer",
        y=y,
        width=0.25,
        hue="model_kmer",
        palette="turbo",
        boxprops={"zorder": 2},
        showfliers=False,
    )
    """for artist in ax.lines:
        artist.set_zorder(10)
    for artist in ax.findobj(PathCollection):
        artist.set_zorder(11)
    ax = sns.stripplot(data=kmer_subset01, x="model_kmer",
        y=y, ax=ax , hue="model_kmer",palette="turbo",jitter=0.4, size=1.5)"""

    if len(selected_kmers) < 6:
        pairs = [
            (selected_kmers[i], selected_kmers[j])
            for i in range(len(selected_kmers))
            for j in range(i + 1, len(selected_kmers))
        ]
    else:
        pairs = [
            (e1, e2)
            for e1, e2 in zip(selected_kmers, selected_kmers[1:] + [selected_kmers[0]])
        ]

    #annotator = Annotator(ax, pairs, data=kmer_subset01, x="Kmers", y="Lengths")
    #annotator.configure(test="Mann-Whitney", text_format="star", loc="inside")
    #annotator.apply_test()
    #ax, test_results = annotator.annotate()

    plt.minorticks_on()
    plt.ylim(0, 45)
    plt.gca().yaxis.grid(True, which="major")
    plt.gca().yaxis.grid(True, which="minor", linewidth=0.2)
    fig.savefig(output_name, dpi=300)


def plot_kmer_noise(df, output_file):
    plot_violin(df, ["A"*9, "C"*9, "G"*9, "T"*9, ], output_file, y="model_stdv")


def plot_kmer_length_ridgeplot(df, output_file):
    df["sample_length"] = df["end_idx"] - df["start_idx"]
    mask = df['sample_length'].values < 40
    df = df.loc[mask]
    plot_boxenplot(df, ["C"*9, "G"*9, "T"*9, "CCAGGCTGG", "TGTGTGTGT", "CCAGCCTGG", "GTGTGTGTG", "TTTTATTTT"], output_file, y="sample_length")
    """df = df.loc[mask]
    selected_kmers = ["A"*9, "C"*9, "G"*9, "T"*9, "CCAGGCTGG", "TGTGTGTGT", "CCAGCCTGG", "GTGTGTGTG"]
    df = df[df.model_kmer.isin(selected_kmers)]

    # Theme
    sns.set_theme(style="white", rc={"axes.facecolor": (0, 0, 0, 0), 'axes.linewidth':2})
    palette = sns.color_palette("Set2", 12)

    # create a grid with a row for each 'Language'
    g = sns.FacetGrid(df, palette=palette, row="model_kmer", hue="model_kmer", aspect=9, height=1.2)

    # map df - Kernel Density Plot of IMDB Score for each Language
    g.map_dataframe(sns.kdeplot, x="sample_length", fill=True, alpha=.5, bw_method=0.5)
    g.map_dataframe(sns.kdeplot, x="sample_length", color='black', bw_method=0.5)

    # function to draw labels
    def label(x, color, label):
        ax = plt.gca() #get current axis
        ax.text(0, .2, label, color='black', fontsize=9,
                ha="left", va="center", transform=ax.transAxes)
    # iterate grid to plot labels
    g.map(label, "model_kmer")

    # adjust subplots to create overlap
    g.fig.subplots_adjust(hspace=-.5)

    # remove subplot titles
    g.set_titles("")

    # remove yticks and set xlabel
    g.set(yticks=[], xlabel="Event length")
    # remove left spine
    g.despine(left=True)
    # set title

    plt.savefig(output_file, dpi=600)"""





def plot_kmer_length(df, output_file, output_file2, output_file3):
    df["sample_length"] = df["end_idx"] - df["start_idx"]
    plot_violin(df, ["A"*9, "C"*9, "G"*9, "T"*9], output_file, y="sample_length")

    plot_violin(df, ["AAAATAAAA", "AAAAATAAA", "AAAAAAATA", "AAAAAAAAT"], output_file2, y="sample_length")

    plot_violin(df, ["GCGCGCGCG", "ATATATATA", "TATATATAT", "CGCGCGCGC"], output_file3, y="sample_length")

def count_missing_9mers(result):
    all_possible_9mers = {"".join(nine_mer): 0 for nine_mer in product("ACGT", repeat=9)}
    missing_9mers = set(all_possible_9mers.keys()) - set(result.keys())
    return len(missing_9mers), len(all_possible_9mers)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Extract and plot frequency distribution of 9-mers from the f5c events file')
    parser.add_argument('input_file', help='Path to the input TSV file')
    parser.add_argument('output_file', help='Path to the output directory')
    args = parser.parse_args()

    # Extract 9-mers from input TSV file
    df, result = extract_9mers(args.input_file, max_rows=100000000) # 100000000

    # Print total number of k-mer sequences
    print("Total number of 9-mers:")
    print(len(result))

    # Print the most common 9-mers and their frequencies
    print("Most common 9-mers:")
    for nine_mer, frequency in Counter(result).most_common(10):
        print(f"{nine_mer}: {frequency}")

    # Count missing 9-mers
    missing_count, possible_count = count_missing_9mers(Counter(result))
    print(f"\nNumber of possible 9-mers: {possible_count}")
    print(f"\nNumber of 9-mers that do not appear at all: {missing_count}")

    outdir = args.output_file
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    # plot k-mer frequency
    plot_frequency_distribution(result, os.path.join(outdir, "kmer-frequency.png"), max_len=100000000)

    # Turn polars to pandas for plotting
    df = df.select(pl.col(["model_kmer","model_stdv", "start_idx", "end_idx"]))
    df_pd = df.to_pandas()

    # Plot relationship k-mer noise
    plot_kmer_noise(df_pd, os.path.join(outdir, "kmer-noise.png"))
    
    # Plot relationship k-mer length
    # plot_kmer_length_ridgeplot(df_pd, os.path.join(outdir, "kmer-length-ridge.png"))


    plot_kmer_length(df_pd, 
                     os.path.join(outdir, "kmer-length.png"), 
                     os.path.join(outdir, "kmer-length-02.png"),
                     os.path.join(outdir, "kmer-length-03.png"),
                    )
