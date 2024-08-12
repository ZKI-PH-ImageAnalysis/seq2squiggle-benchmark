import matplotlib.pyplot as plt
import numpy as np
import matplotlib.colors as mcolors
import scienceplots
import polars as pl
import argparse
import pyslow5

plt.style.use(['nature'])

def plot_main(ax, input_file):
    # Number of different shades of grey
    num_shades = 6
    # Generate shades of grey using Matplotlib
    shades_of_grey = [mcolors.to_rgba(f'lightgrey', alpha=i/num_shades) for i in range(num_shades)]

    patterns = [
        "AAAACTATG",
        "AAACTATGG",
        "AACTATGGT",
        "ACTATGGTA",
        "CTATGGTAG",
        "TATGGTAGA",
        "ATGGTAGAT",
        "TGGTAGATA",
        "GGTAGATAT"
    ]

    results = []
    reader = pl.read_csv_batched(input_file, separator="\t", batch_size=1000000)
    batches = reader.next_batches(100)
    for df in batches:
        # Filter rows where reference_kmer matches any of the patterns
        filtered_df = df.filter(pl.col("model_kmer").is_in(patterns))
        # Append the filtered data to results
        results.append(filtered_df)
    
    # Concatenate all results
    final_df = pl.concat(results)

    # Extract samples and model_mean columns
    final_df = final_df.select(["model_kmer", "samples", "model_mean"])


    # Extract 5 occurrences for each pattern
    selected_rows = []

    for pattern in patterns:
        pattern_df = final_df.filter(pl.col("model_kmer") == pattern).head(4)
        selected_rows.append(pattern_df)

    

    # Concatenate selected rows for all patterns
    final_selection = pl.concat(selected_rows)


    # Extract real data from samples column and group by k-mer
    kmer_signal_map = {pattern: [] for pattern in patterns}
    kmer_model_mean_map = {pattern: [] for pattern in patterns}
    for row in final_selection.iter_rows(named=True):
        kmer_signal_map[row['model_kmer']].append(list(map(float, row['samples'].split(','))))
        kmer_model_mean_map[row['model_kmer']].append(row['model_mean'])

    # Initialize lists to store signals for different currents
    currents = [[] for _ in range(4)]  # Assuming we want four currents

    # Distribute signals among the currents
    for signals in kmer_signal_map.values():
        for i, signal in enumerate(signals):
            currents[i % 4].append(signal)

    # Assign concatenated and padded signals to current variables
    current = currents[0]
    current_01 = currents[1]
    current_02 = currents[2]
    current_03 = currents[3]
    current_flat = [item for sublist in current for item in sublist]
    current_flat_01 = [item for sublist in current_01 for item in sublist]
    current_flat_02 = [item for sublist in current_02 for item in sublist]
    current_flat_03 = [item for sublist in current_03 for item in sublist]
    
    # Simulated time data based on Signal 1 length
    time = np.linspace(0, len(current_flat), len(current_flat))

    events = []
    event_len = []
    for signal, model_mean in zip(current, kmer_model_mean_map.values()):
        # event_level = model_mean
        duration = len(signal)
        event_len.append(duration)
        events.append([model_mean[0]]*duration)
    events_flat = [item for sublist in events for item in sublist]
    
    # Plotting
    ax.plot(time[:len(events_flat)], events_flat, label='Events for Signal 1', color='#d62728')
    ax.plot(time[:len(current_flat)], current_flat, label='Signal 1', color='darkgrey')
    for i, signal in enumerate([current_flat_01, current_flat_02, current_flat_03]):
        ax.plot(time[:len(signal)], signal, label='Signals 2-4' if i == 0 else "", color=shades_of_grey[4 - i])


    # Grey lines
    e_len_sum = 0
    for e_len in event_len:
        e_len_sum += e_len
        ax.axvline(e_len_sum, color='gray', linestyle='--', linewidth=0.5)

    # Single DNA base labels
    labels = ['C', 'T', 'A', 'T', 'G', 'G', 'T', 'A', 'G']
    label_pos = 0 
    for pos, label in zip(event_len, labels):
        start_pos = label_pos
        label_pos += pos
        end_pos = label_pos
        midpoint = (start_pos + end_pos) / 2
        
        ax.text(midpoint, 45, label, ha='center')
    ax.text(-3.5, 45, 'Reference', ha='right')

    # KMer labels
    label_pos = 0 
    for pos, label in zip(event_len, patterns):
        start_pos = label_pos
        label_pos += pos
        end_pos = label_pos
        midpoint = (start_pos + end_pos) / 2
        ax.text(midpoint, 20, label, ha='center', rotation=90, fontsize=5)
    # Add description "Reference" next to the labels
    ax.text(-3.5, 25, 'Kmer', ha='right')

    # Remove x-axis
    ax.set_xlabel('Time')  # Set the x-axis label
    # Hide minor tick labels on the x-axis
    ax.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)

    # Adjusting the plot limits and layout
    ax.set_ylabel('Current (pA)')
    ax.legend(loc='upper left')
    ax.set_ylim(55, 145)  # Adjust y-limits to make space for labels at the bottom

def plot_slow5(ax, input_file, ideal_file):
    # open file
    s5 = pyslow5.Open(input_file,'r')
    currents = []
    # or use directly in a for loop
    for read in s5.seq_reads(pA=True, aux='all'):
        currents.append(read['signal'])
    s5.close()

    # open file
    s6 = pyslow5.Open(ideal_file,'r')
    ideal_currents = []
    # or use directly in a for loop
    for read in s6.seq_reads(pA=True, aux='all'):
        ideal_currents.append(read['signal'])
    s6.close()



    time = np.linspace(0, len(currents[0]), len(currents[0]))
    # Plotting
    shades_of_blue = [mcolors.to_rgba(f'#1f77b4', alpha=i/7) for i in range(1, 8)]

    # ax.plot(time[:len(ideal_currents[0])], ideal_currents[0], label='Predicted Events for Signal 1', color='#d62728')
    # ax.plot(time[:len(currents[0])], currents[0], label='Predicted Signal 1', color=shades_of_blue[6])
    # Plotting Signal 2-4 and summarizing them under one label
    for i in range(0, 4):
        ax.plot(time[:len(currents[i])], currents[i], label='Predicted Signals 1-4' if i == 1 else "", color=shades_of_blue[3 - i])


    ax.set_xlabel('Time')  # Set the x-axis label
    # Hide minor tick labels on the x-axis
    ax.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)

    # Adjusting the plot limits and layout
    ax.set_ylabel('Current (pA)')
    ax.legend(loc='upper left')
    ax.set_ylim(55, 145)  # Adjust y-limits to make space for labels at the bottom

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Process some TSV data and generate a plot.')
    parser.add_argument('input_file', type=str, help='The input TSV file path.')
    parser.add_argument('slow5', type=str, help='The input slow5 file path.')
    parser.add_argument('slow5_ideal', type=str, help='The input slow5 file path withoout noise.')
    parser.add_argument('output_file', type=str, help='The output image file path.')
    args = parser.parse_args()

    fig, axs = plt.subplots(2, 1, figsize=(10, 8))

    plot_main(axs[0], args.input_file)
    plot_slow5(axs[1], args.slow5, args.slow5_ideal)

    axs[0].text(-0.05, 0.99, 'A', transform=axs[0].transAxes, fontsize=12, va='top', ha='left', weight='bold')
    axs[1].text(-0.05, 0.99, 'B', transform=axs[1].transAxes, fontsize=12, va='top', ha='left', weight='bold')

    plt.tight_layout()
    plt.subplots_adjust(hspace=0.5)  
    plt.savefig(args.output_file, dpi=600)