#!/usr/bin/env python3

"""
Plot alignment values from pickles
TODO:
accuracy/identity by alignmentlength
index plot alignment by width.

"""

import os
import pandas as pd
import numpy as np
import pickle
import seaborn as sns
import sys
import humanfriendly

import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.ticker import FuncFormatter
from matplotlib.pylab import savefig
from betaduck.prom_beta_plotter_gen import y_yield_to_human_readable, x_yield_to_human_readable
from pandas.api.types import CategoricalDtype

# Set plot defaults
plot_height = 8.27
plot_width = 11.7
plot_aspect = plot_width / plot_height
bins = 50


def reformat_human_friendly_cov(s):
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
    s += " X"
    return s


def y_cov_to_human_readable(y, position):
    # Convert distribution to base pairs
    if y == 0:
        return 0
    y = round(y, 3)
    s = humanfriendly.format_size(y, binary=False)
    return reformat_human_friendly_cov(s)


def get_pickles(pickle_dir):
    return [pickle_file
            for pickle_file in os.listdir(pickle_dir)
            if pickle_file.endswith(".pickle")]


def get_pickle_data(pickle_files):
    dfs = []
    for pickle_file in pickle_files:
        # ['PAD23566', '6fd51fb7', 'pass', '00004.lambda.pickle']
        flowcell, rand_id, quality, num_id = pickle_file.split("_", 3)

        # Open pickle and generate dataframe
        # Check pickle exists
        if not os.path.isfile(pickle_file):
            print("Warning, cannot find pickle")
            continue
        with open(pickle_file, 'rb') as pickle_h:
            pickle_dict = pickle.load(pickle_h)

        # Add to frame
        df = pd.DataFrame(pickle_dict['read_stats'])

        # add tag
        df['tag'] = pickle_dict['tag']

        # add quality
        df['quality'] = quality

        dfs.append(df)

    # Merge dataframes
    return pd.concat(dfs, axis='rows', sort=False, ignore_index=True)


def plot_dist_split_lambda(df, organism_name, plot_name, attribute='accuracy'):
    # Open up a plotting frame
    g = sns.FacetGrid(df, hue="tag",
                      hue_order=[organism_name, 'lambda'],
                      height=plot_height, aspect=plot_aspect)
    g = g.map(sns.distplot, attribute, bins=bins)

    # Set axis formats
    g.ax.xaxis.set_major_formatter(ticker.PercentFormatter(xmax=1))

    # Generate legend
    g.ax.legend(framealpha=0.5)

    # Set x and y labels
    g.ax.set_title("%s Plot to Human And Lambda" % attribute.capitalize())
    g.ax.set_xlabel("%s (%%)" % attribute)
    g.ax.set_ylabel("")
    g.ax.set_yticks([])

    # Format nicely
    g.fig.tight_layout()

    # Savefig
    savefig("%s.png" % plot_name)


def plot_by_error_type_split_lambda(df, organism_name, plot_name, hue='tag', tag=None):
    alignment_types = ['match', 'mismatch', 'insertion', 'deletion']
    id_vars = ['name', 'tag', 'quality']
    df_melted = pd.melt(df[id_vars + alignment_types],
                        id_vars=id_vars,
                        var_name="alignment_type", value_name="bases")

    # Check tag
    if hue == 'tag' and tag is not None:
        sys.exit("Cant have hue as tag and have only one tag")

    # Set hue order and title / legend names
    if hue == 'tag':
        hue_order = [organism_name, 'lambda']
        plot_title = "Bar Plot of Alignment Type to %s And Lambda" % organism_name.capitalize()
    elif hue == 'quality':
        hue_order = ['pass', 'fail']
        plot_title = "Bar Plot of Alignment Typeto %s And Lambda" % organism_name.capitalize()
    else:
        sys.exit("Unrecognised hue")

    # Generate bar plot
    g = sns.catplot(x="alignment_type",
                    y='bases',
                    hue=hue, hue_order=hue_order,
                    data=df_melted, kind='bar')

    # Set titles, legends
    g.ax.set_title(plot_title)
    g.ax.set_xlabel("Alignment Type")
    g.ax.set_ylabel("Mean bases in Type per Read")

    # Change ylabel to human readable
    g.ax.yaxis.set_major_formatter(FuncFormatter(y_yield_to_human_readable))

    # Add legend
    g.ax.legend([hue.capitalize() for hue in hue_order],
                title=hue.capitalize())

    # Format nicely
    g.fig.tight_layout()

    savefig("%s.png" % plot_name)


def plot_counts_by_chromosome(align_df, plot_name):

    # Initialise plot figure
    fig, ax = plt.subplots(1, figsize=(20, 10))

    # Generate arrays for input
    width = align_df['chrLengthProp']
    height = align_df['cov']
    x_pos = align_df['chrCumPropPoint']
    x_label = align_df['chr'].tolist()

    # Generate bar plot
    ax.bar(x_pos, height, width=width)

    # Reformat y axis
    ax.yaxis.set_major_formatter(FuncFormatter(y_cov_to_human_readable))
    ax.set_ylabel("Coverage", fontsize=14)

    # Reformat x axis
    plt.xticks(x_pos, x_label, rotation=45,
               horizontalalignment='center')
    ax.set_xlabel("Chromosome", fontsize=14)

    ax.set_title("Coverage by chromosome plot", fontsize=20)

    savefig("%s.png" % plot_name)


def plot_alignment_length_by_attribute(df, organism_name, plot_name, attribute='accuracy'):
    """
    Generate a correlation plot of accuracy/identity by alignment length
    :param df:
    :param organism_name:
    :param plot_name:
    :param attribute:
    :return:
    """
    # Reduce to 99th quantile
    max_quantile = 0.99
    max_read_length = df['aln_length'].quantile(max_quantile)
    df = df.query('aln_length < %d' % max_read_length)

    # Seaborn nomenclature for lmplots/regplots are a little different
    sns.set_style('darkgrid')

    g = sns.lmplot(x='aln_length', y=attribute, data=df.query("tag == @organism_name"),
                   x_estimator=np.mean, truncate=True, x_bins=10, scatter_kws={'alpha': 0.1},
                   legend=False)

    # Zero base y-axis
    y_max = df[attribute].mean() * 2
    g.set(ylim=(0, y_max))

    # Set axis labels
    g.set_axis_labels("Alignment Length", "%s (%%)" % attribute.capitalize())

    # Set axis formats
    for ax in g.axes[0]:
        ax.yaxis.set_major_formatter(ticker.PercentFormatter(xmax=1))
        ax.xaxis.set_major_formatter(FuncFormatter(y_yield_to_human_readable))

    # Set title
    g.fig.suptitle("%s over Alignment Length for %s" % (attribute.capitalize(), organism_name))

    # Reduce plot to make room for suptitle
    g.fig.subplots_adjust(top=0.95)

    # Save figure
    savefig("%s.pore_speed.png" % plot_name)
    plt.close('all')


