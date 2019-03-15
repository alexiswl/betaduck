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

import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.ticker import FuncFormatter
from matplotlib.pylab import savefig
from betaduck.prom_beta_plotter_gen import y_yield_to_human_readable


# Set plot defaults
plot_height = 8.27
plot_width = 11.7
plot_aspect = plot_width / plot_height
bins = 50


def get_pickles(pickle_dir):
    return [pickle_file
            for pickle_file in os.listdir(pickle_dir)
            if pickle_file.endswith(".pickle")]


def get_pickle_data(pickle_dir, pickle_files):
    dfs = []
    for pickle_file in pickle_files:
        # ['PAD23566', '6fd51fb7', 'pass', '00004.lambda.pickle']
        flowcell, rand_id, quality, num_id = pickle_file.split("_", 3)

        # Open pickle and generate dataframe
        with open(os.path.join(pickle_dir, pickle_file), 'rb') as pickle_h:
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
    g.ax.set_title("%s Plot to Human And Lambda" % attribute.capitalize)
    g.ax.set_xlabel("Accuracy (%)")
    g.ax.set_ylabel("")
    g.ax.set_yticks([]);

    # Format nicely
    g.fig.tight_layout()

    # Savefig
    savefig("%s.png" % plot_name)


def plot_by_error_type_split_lambda(df, organism_name, plot_name, hue='tag', tag=None):
    alignment_types = ['name', 'match', 'mismatch', 'insertion', 'deletion']
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

