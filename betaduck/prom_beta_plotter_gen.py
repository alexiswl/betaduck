#!/usr/bin/env python3
import os

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use('agg')
from scipy import stats
import matplotlib.pyplot as plt
from matplotlib_venn import venn3
import humanfriendly
from matplotlib.ticker import FuncFormatter
from matplotlib.pylab import savefig

import seaborn as sns
import logging


# Plot yield
def plot_yield(dataset, name, plots_dir):
    """Plot an estimated yield plot and a histogram plot for each sample but by each flowcell"""
    # Plot total yield for the sample
    # Yield plot
    # Set up plotting structure
    fig, ax = plt.subplots(1)

    # Plot setting start_time_float as axis index
    dataset.set_index("start_time_float_by_sample")["yield"].plot(ax=ax)

    # Set x and y ticks
    ax.yaxis.set_major_formatter(FuncFormatter(y_yield_to_human_readable))
    ax.xaxis.set_major_formatter(FuncFormatter(x_yield_to_human_readable))

    # Set x and y labels
    ax.set_title("Yield over time for %s" % name)
    ax.set_xlabel("Time in (HH:MM)")
    ax.set_ylabel("Cumulative Yield")

    # Format nicely
    fig.tight_layout()

    # Save and close figure
    savefig(os.path.join(plots_dir, "%s.yield.png" % name))
    plt.close('all')


# Plot yield by quality
def plot_yield_by_quality(dataset, name, plots_dir):
    # Set up plot
    fig, ax = plt.subplots(1)

    # Iterate through quality and plot each
    q_classes = {"All": "Blue", "Passed": "Green", "Failed": "Red"}
    for quality, col in q_classes.items():
        # Plot the total yield
        if quality == 'All':
            dataset.set_index("start_time_float_by_sample")["yield"].plot(ax=ax, color=col)
        # Plot the yield per quality
        else:
            query = "qualitative_pass == '%s'" % quality
            dataset.set_index("start_time_float_by_sample").query(query)['quality_yield'].plot(ax=ax, color=col)

    # Set x and y ticks
    ax.yaxis.set_major_formatter(FuncFormatter(y_yield_to_human_readable))
    ax.xaxis.set_major_formatter(FuncFormatter(x_yield_to_human_readable))

    # Set x and y labels
    ax.set_title("Yield over time (by quality) for %s" % name)
    ax.set_xlabel("Time in (HH:MM)")
    ax.set_ylabel("Cumulative Yield")

    # Configure legend
    ax.legend(q_classes.keys())

    # Format nicely
    fig.tight_layout()

    # Save and close figure
    savefig(os.path.join(plots_dir, "%s.quality.yield.png" % name))
    plt.close('all')


# Plot reads
def plot_reads(dataset, name, plots_dir):
    """Plot an estimated yield plot and a histogram plot for each sample but by each flowcell"""
    # Plot total number of reads for the sample
    # Set up plotting structure
    fig, ax = plt.subplots(1)

    # Plot setting start_time_float as axis index
    dataset.set_index("start_time_float_by_sample")["read_count"].plot(ax=ax)

    # Set x and y ticks
    ax.yaxis.set_major_formatter(FuncFormatter(y_count_to_human_readable))
    ax.xaxis.set_major_formatter(FuncFormatter(x_yield_to_human_readable))

    # Set x and y labels
    ax.set_title("Read count over time for %s" % name)
    ax.set_xlabel("Time in (HH:MM)")
    ax.set_ylabel("Cumulative read count")

    # Format nicely
    fig.tight_layout()

    # Save and close figure
    savefig(os.path.join(plots_dir, "%s.reads.png" % name))
    plt.close('all')


# Plot read by quality
def plot_read_by_quality(dataset, name, plots_dir):
    # Set up plot
    fig, ax = plt.subplots(1)

    # Iterate through quality and plot each
    q_classes = {"All": "Blue", "Passed": "Green", "Failed": "Red"}
    for quality, col in q_classes.items():
        # Plot the total yield
        if quality == 'All':
            dataset.set_index("start_time_float_by_sample")["read_count"].plot(ax=ax, color=col)
        # Plot the yield per quality
        else:
            query = "qualitative_pass == '%s'" % quality
            dataset.set_index("start_time_float_by_sample").query(query)['quality_count'].plot(ax=ax, color=col)

    # Set x and y ticks
    ax.yaxis.set_major_formatter(FuncFormatter(y_count_to_human_readable))
    ax.xaxis.set_major_formatter(FuncFormatter(x_yield_to_human_readable))

    # Set x and y labels
    ax.set_title("Read count over time (by quality) for %s" % name)
    ax.set_xlabel("Time in (HH:MM)")
    ax.set_ylabel("Cumulative Read Count")

    # Configure legend
    ax.legend(q_classes.keys())

    # Format nicely
    fig.tight_layout()

    # Save and close figure
    savefig(os.path.join(plots_dir, "%s.quality.reads.png" % name))
    plt.close('all')


def plot_flowcell(dataset, name, plots_dir):
    # Set up plots 
    fig, ax = plt.subplots()
    fig.set_size_inches(15, 7)

    # Use the formatter we used for the yield plots.
    formatter_y = FuncFormatter(y_yield_to_human_readable)

    # Create the values that make up the numbers on the far-right column of the grid.
    c_w = 10
    c_l = 25
    c_num = 12

    # Create the array
    channels_by_order_array = np.array([[c_no * c_w * c_l + c_w * l_no + w_no + 1
                                         for c_no in np.arange(c_num)
                                         for w_no in np.arange(c_w)]
                                        for l_no in np.arange(c_l)])

    # Use the minknow_column_order function which reference the far-right column for a given row

    # to fill in the rest of the values for each row.
    channels_by_yield_array = np.zeros(channels_by_order_array.shape)

    # Sum the values for each channel.
    channels_by_yield_df = pd.DataFrame(dataset.groupby("channel")['channel_yield'].max())

    # Reset the index and have channel as a column instead of the index.
    channels_by_yield_df.reset_index(level='channel', inplace=True)

    # Iterate through each row of the yield by channel dataframe.
    for yield_row in channels_by_yield_df.itertuples():
        channel_index = [(ix, iy)
                         for ix, row in enumerate(channels_by_order_array)
                         for iy, i in enumerate(row)
                         if int(i) == int(yield_row.channel)][0]

        # Assign channel yield to position in MinKNOW
        channels_by_yield_array[channel_index] = yield_row.channel_yield

    # Plot heatmap
    sns.heatmap(channels_by_yield_array,
                # Remove labels from side, they're not useful in this context.
                xticklabels=False,
                yticklabels=False,
                ax=ax,
                # Prevent extreme values from over-scaling the sidebar.
                robust=True,
                # Use the greens scale but in reverse, similar to MinKNOW.
                cmap=sns.diverging_palette(210, 120, l=55, as_cmap=True),
                # Format keyword args for the side bar.
                cbar_kws={"format": formatter_y,
                          "label": "Bases per channel"})

    # Create three lines down the middle as shown in PromethION MinKNOW.
    [ax.axvline([x], color='white', lw=5) for x in [30, 60, 90]]

    # Nice big title!
    ax.set_title("Map of Yield by Channel for %s" % name, fontsize=25)

    # Ensure labels are not missed.
    fig.tight_layout()

    # Save and close figure
    savefig(os.path.join(plots_dir, "%s.flowcellmap.png" % name))
    plt.close('all')


# Plot histogram
def plot_weighted_hist(dataset, name, plots_dir):
    # Set globals
    num_bins = 50

    # Open up plotting frame
    fig, ax = plt.subplots(1)

    # Get linspacing of histogram
    bins = np.linspace(start=0, stop=dataset['sequence_length_template'].max(), num=num_bins)

    # Develop the bin width
    bin_width = bins[1] - bins[0]

    # Plot weighted histogram
    dataset['sequence_length_template'].plot(kind="hist", ax=ax, density=1, bins=bins,
                                             alpha=0.6, weights=dataset['sequence_length_template'])

    # Set the axis formatters
    def y_hist_to_human_readable_seq(y, position):
        # Convert distribution to base pairs
        if y == 0:
            return 0
        s = humanfriendly.format_size(bin_width * dataset['sequence_length_template'].sum() * y, binary=False)
        return reformat_human_friendly(s)

    ax.yaxis.set_major_formatter(FuncFormatter(y_hist_to_human_readable_seq))
    ax.xaxis.set_major_formatter(FuncFormatter(x_hist_to_human_readable))

    # Set titles
    ax.set_title("Read Distribution Graph for %s" % name)
    ax.grid(color='black', linestyle=':', linewidth=0.5)
    ax.set_ylabel("Bases per bin")
    ax.set_xlabel("Read length")

    # Ensure labels are not missed.
    fig.tight_layout()

    # Save and close figure
    savefig(os.path.join(plots_dir, "%s.weighted.hist.png" % name))
    plt.close('all')


def plot_read_hist(dataset, name, plots_dir):
    # Much simpler histogram with seaborn
    # Set max quantile, we need to reduce this for the read length histogram as it's not weighter
    max_quantile = 0.99
    max_read_length = dataset['sequence_length_template'].quantile(max_quantile)
    trimmed = dataset.query('sequence_length_template < %d' % max_read_length)
    # Open up a plotting frame
    fig, ax = plt.subplots(1)

    # Set seaborn style
    sns.set_style("darkgrid")

    # Plot distribution
    sns.distplot(trimmed['sequence_length_template'],
                 hist=True, kde=True, ax=ax)

    # Despine left axis
    sns.despine(fig=fig, ax=ax, left=True)

    # Set titles
    ax.set_title("Read Length Distribution")

    # Set x and y lables
    ax.set_xlabel("Read Length")

    # Set x-axis
    ax.xaxis.set_major_formatter(FuncFormatter(x_hist_to_human_readable))

    # Ensure labels are not missed
    fig.tight_layout()

    # Save and close figure
    savefig(os.path.join(plots_dir, "%s.unweighted.hist.png" % name))


def plot_quality_hist(dataset, name, plots_dir):
    # Much simpler histogram with seaborn

    # Open up a plotting frame
    fig, ax = plt.subplots(1)

    # Set seaborn style
    sns.set_style("darkgrid")

    # Plot distribution
    sns.distplot(dataset['mean_qscore_template'],
                 hist=True, kde=True, ax=ax)

    # Despine left axis
    sns.despine(fig=fig, ax=ax, left=True)

    # Set titles
    ax.set_title("Mean QScore Distribution")

    # Set x and y lables
    ax.set_xlabel("Mean QScore")

    # Ensure labels are not missed
    fig.tight_layout()

    # Save and close figure
    savefig(os.path.join(plots_dir, "%s.quality.hist.png" % name))


def reformat_human_friendly(s):
    """
    humanfriendly module returns with a few quirks
    1 = 1 byte ==> 1 b
    2 = 2 bytes ==> 2 bytes
    1000 = 1 KB ==> 1 Kb
    """
    s = s.replace(" byte", "")
    s = s.replace(" bytes", "")
    s = s.replace("B", "")
    s = s.replace("s", "")
    s += "b"
    return s


def y_yield_to_human_readable(y, position):
    # Convert distribution to base pairs
    if y == 0:
        return 0
    y = round(y, 3)
    s = humanfriendly.format_size(y, binary=False)
    return reformat_human_friendly(s)


def y_count_to_human_readable(y, position):
    # Use the same as y yield but strip off the last 'b'
    if y == 0:
        return 0
    y = round(y, 3)
    s = humanfriendly.format_size(y, binary=False)
    s = reformat_human_friendly(s).rstrip('b')
    return s


def x_yield_to_human_readable(x, position):
    # Convert time in seconds to hours or minutes
    hours = int(x // 3600)
    minutes = int((x % 3600) // 60)
    seconds = int(x % 60)
    if x == 0:
        return 0
    s = f"{hours:02d}:{minutes:02d}"
    return s


def x_hist_to_human_readable(x, position):
    # Convert distribution to base pairs
    if x == 0:
        return 0
    s = humanfriendly.format_size(x, binary=False)
    return reformat_human_friendly(s)


def plot_events_ratio(dataset, name, plots_dir):
    # Seaborn nomenclature for reg/lm plots are a little different

    # Set the background style for the plot
    sns.set_style('darkgrid')

    # Generate the plot 
    g = sns.lmplot(x='start_time_float_by_sample', y='events_ratio', data=dataset,
                   hue='qualitative_pass', hue_order=['Passed', 'Failed'],
                   x_estimator=np.mean, truncate=True, x_bins=10, scatter_kws={'alpha': 0.1},
                   legend=False)

    # Create legend and rename
    leg_title = "Read Quality"
    leg = g.ax.legend(title=leg_title, framealpha=0.5)
    for lh in leg.legendHandles:
        lh.set_alpha(1)

    # Zero base y-axis
    y_max = dataset['events_ratio'].mean() * 2
    g.set(ylim=(0, y_max))

    # Set x and y labels
    g.set_axis_labels("Time in (HH:MM)", "Events ratio (events / base)")

    # Set title
    g.fig.suptitle("Events Ratio Graph for %s" % name)

    # Set x and y ticks:
    for ax in g.axes[0]:
        ax.xaxis.set_major_formatter(FuncFormatter(x_yield_to_human_readable))
        ax.yaxis.set_major_formatter(FuncFormatter(y_count_to_human_readable))

    # Format nicely
    g.fig.tight_layout()

    # Reduce the plot size to make way for the suptitle
    g.fig.subplots_adjust(top=0.95)

    # Save and close figure
    savefig(os.path.join(plots_dir, "%s.events_ratio.png" % name))
    plt.close("all")


def plot_quality_per_speed(dataset, name, plots_dir):
    # Seaborn nomenclature for joint plots are a little different
    sns.set_style("dark")
    g = sns.jointplot(x='pore_speed', y='mean_qscore_template',
                      data=dataset, kind='hex')

    # Add pearson stat
    g.annotate(stats.pearsonr)

    # Set axis labels
    g.set_axis_labels("Pore Speed (b/s)", "Mean Q-score")

    # Set title
    g.fig.suptitle("Pore Speed vs Q-score for %s" % name)

    # Format nicely.
    g.fig.tight_layout()

    # Reduce plot to make room for suptitle
    g.fig.subplots_adjust(top=0.95)

    # Save and close the figure
    savefig(os.path.join(plots_dir, "%s.speed_vs_qscore.png" % name))
    plt.close('all')


def plot_venn_diagram_of_filtered_data(dataset, filter_dict, name, plots_dir):
    # Iterate through each of the different objects obtaining the shape required for the venn diagram
    time_subset = dataset.query(' & '.join([filter_dict['Time'][0], filter_dict['Events Ratio'][1], filter_dict['Max Read Length'][1]])).shape[0]
    events_subset = dataset.query(' & '.join([filter_dict['Time'][1], filter_dict['Events Ratio'][0], filter_dict['Max Read Length'][1]])).shape[0]
    time_and_events_subset = dataset.query(' & '.join([filter_dict['Time'][1], filter_dict['Events Ratio'][1], filter_dict['Max Read Length'][0]])).shape[0]
    length_subset = dataset.query(' & '.join([filter_dict['Time'][1], filter_dict['Events Ratio'][1], filter_dict['Max Read Length'][0]])).shape[0]
    time_and_length_subset = dataset.query(' & '.join([filter_dict['Time'][0], filter_dict['Events Ratio'][1], filter_dict['Max Read Length'][0]])).shape[0]
    events_and_length_subset = dataset.query(' & '.join([filter_dict['Time'][1], filter_dict['Events Ratio'][0], filter_dict['Max Read Length'][0]])).shape[0]
    all_subset = dataset.query(' & '.join([filter_dict['Time'][0], filter_dict['Events Ratio'][0], filter_dict['Max Read Length'][0]])).shape[0]
    fig, ax = plt.subplots()
    venn3(subsets=(time_subset, events_subset, time_and_events_subset, length_subset,
                   time_and_length_subset, events_and_length_subset, all_subset),
          set_labels=['Time', "Events Ratio", "Max Read Length"],
          ax=ax)

    # Set titles
    ax.set_title("Reads excluded by condition")

    # Ensure labels are not missed
    fig.tight_layout()

    # Save and close figure
    savefig(os.path.join(plots_dir, "%s.venn_diagram.png" % name))

def plot_quality_per_readlength(dataset, name, plots_dir):
    # Set max quantile, we need to reduce this for the read length histogram as it's not weighter
    max_quantile = 0.98
    max_read_length = dataset['sequence_length_template'].quantile(max_quantile)
    trimmed = dataset.query('sequence_length_template < %d' % max_read_length)

    # Seaborn nomenclature for joint plots are a little different
    sns.set_style("dark")

    g = sns.jointplot(x='sequence_length_template', y='mean_qscore_template',
                      data=trimmed, kind='hex')

    # Add pearson stat
    g.annotate(stats.pearsonr)

    # Set axis labels
    g.set_axis_labels("Sequence length", "Mean Q-score")

    # Set x ticks: 
    g.ax_joint.xaxis.set_major_formatter(FuncFormatter(x_hist_to_human_readable))

    # Set title
    g.fig.suptitle("Sequence length against Q-score for %s" % name)

    # Format nicely.
    g.fig.tight_layout()

    # Reduce plot to make room for suptitle
    g.fig.subplots_adjust(top=0.95)

    # Save and close the figure
    savefig(os.path.join(plots_dir, "%s.length_vs_qscore.png" % name))
    plt.close('all')


def plot_pair_plot(dataset, name, plots_dir):
    # Globals used
    sample_size = 10000

    # Plot everything side by side.
    sns.set_style("darkgrid")

    # Select columns to plot
    items = ["mean_qscore_template", "pore_speed",
             "sequence_length_template", "events_ratio",
             "qualitative_pass"]

    # Changing the axis labels is easier done first.
    rename_columns = {"mean_qscore_template": "Mean QScore Template",
                      "pore_speed": "Pore Speed (b/s)",
                      "sequence_length_template": "Read Length",
                      "events_ratio": "Events / base"}

    # We need to further lower the events ratio
    events_ratio_threshold = 5

    # Get sample
    sample_set = dataset.query("events_ratio < %d" % events_ratio_threshold).filter(items=items).sample(sample_size)

    # Plot grid
    g = sns.PairGrid(sample_set.rename(columns=rename_columns))

    # Scatter in the top corner
    g.map_upper(plt.scatter, s=1, alpha=0.5)  # hue='qualitative_pass', hue_order=["Passed", "Failed"])

    # Kde plots in the bottom corner
    g.map_lower(sns.kdeplot, shade=True, shade_lowest=False)

    # Distribution plots down the middle
    g.map_diag(plt.hist, weights=sample_set['sequence_length_template'])

    # Use the funcformmatter for the Read Lengths (3rd row and column (so 2 on a 0 indexed array))
    g.axes[2, 2].xaxis.set_major_formatter(FuncFormatter(y_yield_to_human_readable))
    g.axes[2, 2].yaxis.set_major_formatter(FuncFormatter(y_yield_to_human_readable))

    # Set title
    g.fig.suptitle("Pair plot for %s" % name)

    # Set FuncFormatter on ax_joint (later).

    # Reduce plot to make room for suptitle
    g.fig.subplots_adjust(top=0.95)

    # Save figure
    savefig(os.path.join(plots_dir, "%s.pair_plot.png" % name))
    plt.close('all')


def plot_quality_over_time(dataset, name, plots_dir):
    # Set global number of bins
    bins = 16

    # Set max and min values for times
    min_value = 0
    max_value = max(dataset['start_time_float_by_sample'])
    # Generate the cut
    dataset['start_time_float_by_sample_bin'] = pd.cut(dataset['start_time_float_by_sample'], bins=bins)
    dataset['start_time_float_by_sample_bin_left'] = dataset['start_time_float_by_sample_bin'].apply(lambda x: max(x.left, min_value))
    dataset['start_time_float_by_sample_bin_right'] = dataset['start_time_float_by_sample_bin'].apply(lambda x: min(x.right, max_value))
    dataset['start_time_float_by_sample_bin_str'] = dataset.apply(lambda x: ' - '.join(map(str, [x_yield_to_human_readable(x.start_time_float_by_sample_bin_left, None),
                                                                                                 x_yield_to_human_readable(x.start_time_float_by_sample_bin_right, None)])),
                                                                  axis='columns')

    # Generate a ridges plot, splitting the dataframe into fifteen bins.
    # Initialize the FacetGrid object
    pal = sns.cubehelix_palette(bins, rot=-.25, light=.7)
    g = sns.FacetGrid(dataset,
                      row="start_time_float_by_sample_bin_str", hue="start_time_float_by_sample_bin_str",
                      aspect=15, height=2, palette=pal)

    # Draw the densities in a few steps
    g.map(sns.kdeplot, "mean_qscore_template", clip_on=False, shade=True, alpha=1, lw=1.5, bw='scott')
    g.map(sns.kdeplot, "mean_qscore_template", clip_on=False, color="w", lw=2, bw='scott')
    g.map(plt.axhline, y=0, lw=2, clip_on=False)
    g.map(plt.axhline, y=0, lw=2, clip_on=False)

    # Define and use a simple function to label the plot in axes coordinates
    # Not sure why, but input must be name as function.
    def set_xlabel(color, label):
        ax = plt.gca()
        ax.text(0, .6, label, fontweight="bold", color=color,
                ha="left", va="center", transform=ax.transAxes)

    g.map(set_xlabel, color='black', label="start_time_float_by_sample_bin_str")

    # Remove axes details that don't play well with overlap
    g.set_titles("")
    g.set(yticks=[])
    g.despine(bottom=True, left=True)

    # Add title
    g.fig.suptitle("Quality distribution over time")
    # Reset xlabel
    g.set_xlabels("Mean Q-Score")

    # Reduce plot to make room for suptitle
    g.fig.subplots_adjust(top=0.95)

    # Save figure
    savefig(os.path.join(plots_dir, "%s.q_score.time.split.png" % name))
    plt.close('all')


def plot_pore_speed(dataset, name, plots_dir):
    # Plot setting start_time_float as axis index

    # Seaborn nomenclature for lmplots/regplots are a little different
    sns.set_style('darkgrid')

    g = sns.lmplot(x='start_time_float_by_sample', y='pore_speed', data=dataset,
                   hue='qualitative_pass', hue_order=["Passed", "Failed"],
                   x_estimator=np.mean, truncate=True, x_bins=10, scatter_kws={'alpha': 0.1},
                   legend=False)

    # Create legend with new alpha
    leg_title = "Read Quality"
    leg = g.ax.legend(title=leg_title, framealpha=0.5)
    for lh in leg.legendHandles:
        lh.set_alpha(1)

    # Zero base y-axis
    y_max = dataset['pore_speed'].mean() * 2
    g.set(ylim=(0, y_max))

    # Set axis labels
    g.set_axis_labels("Time (HH:MM)", "Pore Speed (bases / second)")

    # Set axis formats
    for ax in g.axes[0]:
        ax.xaxis.set_major_formatter(FuncFormatter(x_yield_to_human_readable))

    # Set title
    g.fig.suptitle("Pore speed over time")

    # Reduce plot to make room for suptitle
    g.fig.subplots_adjust(top=0.95)

    # Save figure
    savefig(os.path.join(plots_dir, "%s.pore_speed.png" % name))
    plt.close('all')


def print_stats(dataset, name, plots_dir):
    percentiles = [0.1, 0.25, 0.5, 0.75, 0.9]
    # Get total yield
    total_bp = dataset['sequence_length_template'].sum()
    total_bp_h = reformat_human_friendly(humanfriendly.format_size(total_bp, binary=False))
    # Describe length
    length_describe = dataset['sequence_length_template'].describe(percentiles=percentiles).to_string()
    # Describe quality
    qual_describe = dataset['mean_qscore_template'].describe(percentiles=percentiles).to_string()

    # Reformat each of the methods such that they're rounded to two decimal places
    length_describe = '\n'.join([qual_line.split()[0].ljust(8) + "\t" +
                                 "{:21.2f}".format(float(qual_line.split()[1]))
                                 for qual_line in length_describe.split("\n")])
    qual_describe = '\n'.join([qual_line.split()[0].ljust(8) + "\t" +
                               "{:21.2f}".format(float(qual_line.split()[1]))
                               for qual_line in qual_describe.split("\n")])

    # Calculate the N50:
    nx = []
    seq_length_sorted_as_series = dataset['sequence_length_template'].sort_values().reset_index(drop=True)
    seq_length_cumsum_as_series = seq_length_sorted_as_series.cumsum()
    for index, seq_value in seq_length_sorted_as_series.iteritems():
        if (seq_length_cumsum_as_series[index] <= total_bp * percentiles[len(nx)] <=
                seq_length_cumsum_as_series[index + 1]):
            nx.append(seq_value)
        if len(nx) == len(percentiles):
            # Found all the percentiles, no need to continue.
            break
    nx_h = [reformat_human_friendly(humanfriendly.format_size(n_x_value, binary=False))
            for n_x_value in nx]

    # Get run duration, from first read to last read.
    duration = dataset["start_time_float_by_sample"].max()  # In seconds
    hours, remainder = divmod(duration, 3600)
    minutes, seconds = divmod(remainder, 60)
    run_duration_h = f"{hours} hours, {minutes} minutes, {seconds:2,.0f} seconds"

    # Print these stats
    sample_name = dataset["sample_id"].unique().item()
    with open(os.path.join(plots_dir, "%s.stats.txt" % name), 'a') as output_handle:
        # Print total basepairs
        output_handle.write("# Stats for sample '%s' #\n" % sample_name)
        output_handle.write("Total basepairs:\n")
        output_handle.write(f"\t{total_bp:16,.0f}\t|\t{total_bp_h.rjust(9)}\n")
        output_handle.write("Description of Read Lengths:\n")
        # Tab indent each of the descriptor lines
        output_handle.writelines(f"\t{len_line}\n"
                                 for len_line in length_describe.split("\n"))
        output_handle.write("Description of Read Qualities:\n")
        # Tab indent each of the descriptor lines
        output_handle.writelines(f"\t{qual_line}\n"
                                 for qual_line in qual_describe.split("\n"))
        output_handle.write("NX values:\n")
        output_handle.writelines(f"\tN{100*percentile:02.0f}:\t{nx_value:8,.0f}\t|\t{nx_h_value.rjust(9)}\n"
                                 for percentile, nx_value, nx_h_value in zip(percentiles, nx, nx_h))
        output_handle.write(f"\t{duration:8,.1f} seconds\t|\t{run_duration_h}\n")


def plot_data(dataset, name, plots_dir):
    # Plot things
    # Matplotlib base plots
    plotting_functions = [plot_yield, plot_yield_by_quality, plot_reads, plot_read_by_quality,
                          plot_weighted_hist, plot_read_hist, plot_flowcell, plot_pore_speed,
                          plot_quality_hist, plot_quality_over_time,
                          plot_quality_per_speed, plot_quality_per_readlength,
                          plot_events_ratio, plot_pair_plot]

    # Just iterate through each of the plotting methods.
    for function in plotting_functions:
        function(dataset=dataset, name=name, plots_dir=plots_dir)

    # Print out stats
    logging.info("Finishing plotting")
