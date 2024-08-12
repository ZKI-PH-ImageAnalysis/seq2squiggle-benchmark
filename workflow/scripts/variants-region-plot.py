import argparse
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

def read_data(file_path):
    """ Function to read data from the input text file """
    data = {
        'False positives': {
            'homopolymers': 0,
            'STRs': 0,
            'other regions': 0
        },
        'False negatives': {
            'homopolymers': 0,
            'STRs': 0,
            'other regions': 0
        },
        'True positives': {
            'homopolymers': 0,
            'STRs': 0,
            'other regions': 0
        }
    }

    with open(file_path, 'r') as f:
        for line in f:
            line = line.strip()
            if line:
                category, value = line.split(': ')
                if 'False positives' in category:
                    group, subtype = category.split(' due to ')
                    data['False positives'][subtype] = int(value)
                elif 'False negatives' in category:
                    group, subtype = category.split(' due to ')
                    data['False negatives'][subtype] = int(value)
                elif 'True positives' in category:
                    group, subtype = category.split(' due to ')
                    data['True positives'][subtype] = int(value)

    return data

def plot_grouped_barplot(tool1_data, tool2_data, output_file):
    """ Function to plot the grouped barplot and save as PNG """
    sns.set_theme(style="whitegrid")
    sns.set_context("paper", font_scale=1.5)
    
    categories = ['homopolymers', 'STRs', 'other regions']
    labels = ['False positives', 'False negatives', 'True positives']
    x = np.arange(len(labels))  # the label locations
    width = 0.35  # the width of the bars

    fig, ax = plt.subplots(figsize=(10, 6))

    # Extract data for tool1 and tool2
    tool1_counts = [tool1_data[label][category] for label in labels for category in categories]
    tool2_counts = [tool2_data[label][category] for label in labels for category in categories]


    # Reshape tool1_counts and tool2_counts to match the shape needed for plotting
    tool1_counts = [tool1_counts[i:i+len(categories)] for i in range(0, len(tool1_counts), len(categories))]
    tool2_counts = [tool2_counts[i:i+len(categories)] for i in range(0, len(tool2_counts), len(categories))]

    # Plotting with customized colors
    colors_seq2squiggle = sns.color_palette('Blues', len(labels))
    colors_squigulator = sns.color_palette('Oranges', len(labels))

    rects1 = []
    rects2 = []
    for i in range(len(labels)):
        rects1.append(ax.bar(x - width/2, tool1_counts[i], width, bottom=np.sum(tool1_counts[:i], axis=0), label=labels[i], color=colors_seq2squiggle[i]))
        rects2.append(ax.bar(x + width/2, tool2_counts[i], width, bottom=np.sum(tool2_counts[:i], axis=0), label=labels[i], color=colors_squigulator[i]))


    # Add some text for labels, title and custom x-axis tick labels, etc.
    ax.set_ylabel('Counts')
    ax.set_title('Comparison of seq2squiggle and squigulator')
    ax.set_xticks(x)
    ax.set_xticklabels(categories)
    ax.legend(title='Categories', loc='upper right')

    fig.tight_layout()

    # Save the plot as PNG file
    plt.savefig(output_file, dpi=600)
    plt.close()


def plot_stacked_barplot(tool1_data, tool2_data, output_file):
    """ Function to plot the stacked barplot and save as PNG """
    sns.set_theme(style="whitegrid")
    sns.set_context("paper", font_scale=1.5)
    
    categories = ['homopolymers', 'STRs', 'other regions']
    labels = ['False positives', 'False negatives', 'True positives']
    x = np.arange(len(categories))  # the label locations
    width = 0.35  # the width of the bars

    fig, ax = plt.subplots(figsize=(14, 7))

    # Extract data for tool1 and tool2 and convert to percentages
    tool1_counts = [[tool1_data[label][category] for category in categories] for label in labels]
    tool2_counts = [[tool2_data[label][category] for category in categories] for label in labels]

    tool1_sums = [sum(category) for category in zip(*tool1_counts)]
    tool2_sums = [sum(category) for category in zip(*tool2_counts)]

    tool1_percents = [[(count / total) * 100 if total != 0 else 0 for count, total in zip(counts, tool1_sums)] for counts in tool1_counts]
    tool2_percents = [[(count / total) * 100 if total != 0 else 0 for count, total in zip(counts, tool2_sums)] for counts in tool2_counts]

    # Plotting with customized colors
    colors_seq2squiggle = sns.color_palette('Blues', len(labels))
    colors_squigulator = sns.color_palette('Oranges', len(labels))

    bottoms_tool1 = np.zeros(len(categories))
    bottoms_tool2 = np.zeros(len(categories))

    for i in range(len(labels)):
        ax.bar(x - width/2, tool1_percents[i], width, bottom=bottoms_tool1, label=f'seq2squiggle {labels[i]}', color=colors_seq2squiggle[i])
        ax.bar(x + width/2, tool2_percents[i], width, bottom=bottoms_tool2, label=f'squigulator {labels[i]}', color=colors_squigulator[i])
        bottoms_tool1 += np.array(tool1_percents[i])
        bottoms_tool2 += np.array(tool2_percents[i])

    # Add some text for labels, title and custom x-axis tick labels, etc.
    ax.set_ylabel('Percentage')
    # ax.set_title('Comparison of seq2squiggle and squigulator')
    ax.set_xticks(x)
    ax.set_xticklabels(categories)
    
    # Move the legend below the plot
    ax.legend(title='Categories', loc='upper center', bbox_to_anchor=(0.5, -0.15), fontsize='small', ncol=3)

    fig.tight_layout()

    # Save the plot as PNG file
    plt.savefig(output_file, dpi=600)
    plt.close()


def main():
    # Argument parsing
    parser = argparse.ArgumentParser(description='Generate grouped barplot comparing tool results.')
    parser.add_argument('tool1_file', type=str, help='Path to seq2squiggle results text file')
    parser.add_argument('tool2_file', type=str, help='Path to squigulator results text file')
    parser.add_argument('output_file', type=str, help='Output PNG file for the barplot')
    args = parser.parse_args()

    # Read data from the input files
    tool1_data = read_data(args.tool1_file)
    tool2_data = read_data(args.tool2_file)
    
    # Plot and save the grouped barplot as PNG
    plot_stacked_barplot(tool1_data, tool2_data, args.output_file)
    print(f"Stacked barplot saved as '{args.output_file}'")

if __name__ == "__main__":
    main()