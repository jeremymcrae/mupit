"""
Copyright (c) 2016 Genome Research Ltd.

Permission is hereby granted, free of charge, to any person obtaining a copy of
this software and associated documentation files (the "Software"), to deal in
the Software without restriction, including without limitation the rights to
use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies
of the Software, and to permit persons to whom the Software is furnished to do
so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
"""

from pandas import Series
from numpy import log10

import matplotlib
matplotlib.use("Agg")
from matplotlib.backends.backend_pdf import PdfPages

from matplotlib import pyplot

def plot_enrichment(enriched, num_tests, path, chrom="chrom", symbol="hgnc",
    position="min_pos", p_columns=None):
    """ make Manhattan plots for loss-of-function and functional variants
    
    Args:
        enriched: data frame containing columns for chr, coord (position), and
             P-values from testing for enrichment with loss-of-function and
             functional consequence variants.
        num_tests: number of tests performed (used by Bonferroni and FDR
             correction).
        path: path to save the plot to.
        chrom: column name for chromosome labels
        symbol: column name for gene symbols
        position: column name for nucleotide position coordinate
        p_columns: list of columns of p-values to plot (one page of pdf per
            p-value column)
    """
    
    if p_columns is None:
        raise ValueError("No p-value columns defined!")
    
    if type(p_columns) == str:
        p_columns = [p_columns]
    
    enriched = enriched[[symbol, chrom, position] + p_columns].copy()
    
    # recode the sex chromosomes as integers
    enriched[chrom] = enriched[chrom].astype(str)
    enriched[chrom][enriched[chrom] == "X"] = "23"
    enriched[chrom] = enriched[chrom].astype(int)
    
    # drop rows with missing values
    for column in enriched.columns:
        enriched = enriched[~enriched[column].isnull()]
    
    # sort the table by genome coordinates
    enriched = enriched.reset_index()
    enriched = enriched.sort([chrom, position])
    
    enriched["pos"] = adjust_coordinates(enriched, chrom, position)
    
    with PdfPages(path) as pdf:
        for column in p_columns:
            make_figure(enriched, num_tests, column, pdf, chrom)

def make_figure(enriched, num_tests, test, pdf, chrom="chrom", symbol="hgnc"):
    """
    
    Args:
        enriched: data frame containing columns for chr, coord (position), and
             P-values from testing for enrichment with loss-of-function and
             functional consequence variants.
        num_tests: number of tests performed (used by Bonferroni and FDR
             correction).
        test: column name for p-values for the test to be plotted
        pdf: PdfPages plotting device, for optional multipage output
        chrom: column name for chromosome labels
        symbol: column name for gene symbols
    """
    
    fig = pyplot.figure(figsize=(10, 6))
    title = fig.suptitle(test)
    ax = fig.gca()
    
    ticks = get_ticks(enriched, chrom)
    enriched["log10_p"] = -log10(enriched[test])
    threshold = 0.01/num_tests
    
    # plot each chromosome one by one, in alternating colors
    colors = ["lightblue", "darkblue"]
    color = colors[0]
    for key in sorted(set(enriched[chrom])):
        points = enriched[enriched[chrom] == key]
        line2d = pyplot.plot(points["pos"], points["log10_p"],
            axes=ax, linestyle='None', marker="o", markeredgewidth=0,
            markeredgecolor=color, markerfacecolor=color)
        color = list(set(colors) - set([color]))[0]
    
    # expand the axis limits slightly
    xlim = ax.set_xlim([min(enriched["pos"]) * -0.03, max(enriched["pos"]) * 1.03])
    ylim = ax.set_ylim([0, max(enriched["log10_p"]) * 1.03])
    
    # add a horizontal line to indicate genomewide significance, and annotate
    # all the genes which exceed that threshold
    hline = ax.axhline(y=-log10(threshold), linestyle="dashed", color="red")
    signif = enriched[enriched["log10_p"] > -log10(threshold)]
    for (key, row) in signif.iterrows():
        text = ax.text(row["pos"], row["log10_p"], row[symbol],
            fontsize="xx-small", fontstyle="italic")
    
    annotate_plot(ax, ticks, chrom)
    
    pdf.savefig()

def annotate_plot(ax, ticks, chrom="chrom"):
    """ annotate a manhattan plot with axis labels, tickmarks etc
    
    Args:
        ax: matplotlib Axes for the current plot
        chrom: column name for chromosome labels
        ticks: list of positions for chromosome tick marks.
    """
    
    # Hide the right, top and bottom spines
    e = ax.spines['right'].set_visible(False)
    e = ax.spines['top'].set_visible(False)
    e = ax.spines['bottom'].set_visible(False)
    e = ax.tick_params(direction='out')
    
    # Only show ticks on the left and bottom spines
    e = ax.yaxis.set_ticks_position('left')
    e = ax.xaxis.set_ticks_position('bottom')
    
    # add in the chromosome ticks and labels
    e = ax.set_xticks(ticks)
    e = ax.set_xticklabels(sorted(set(enriched[chrom])))
    
    # define the axes labels
    e = ax.set_xlabel("Chromosome")
    e = ax.set_ylabel("-log10(P)")

def adjust_coordinates(coords, chrom="chrom", position="min_pos"):
    """ get sequential positions across successive chromosomes
    
    Manhattan plots have successive chromosomes plotted one after another. But
    nucleotide coordinates reset between chromosomes. We need to create a new
    coordinate, where the position takes into account the chromosomes that have
    come before.
    
    Args:
        coords: pandas DataFrame containing rows with genome coordinates
        chrom: column name for chromosome labels
        position: column name for nucleotide position coordinates
    
    Returns:
        pandas Series of adjusted positions.
    """
    
    pos = Series([0] * len(coords))
    
    for key in sorted(set(coords[chrom])):
        rows = coords[chrom] == key
        pos[rows] = coords[position][rows] + max(pos)
    
    return pos

def get_ticks(coords, chrom="chrom"):
    """ identify the positions for tick marks for chromosomes
    
    Args:
        coords: pandas DataFrame containing rows with genome coordinates
        chrom: column name for chromosome labels
    
    Returns:
        list of positions for chromosome tick marks.
    """
    
    ticks = []
    for key in sorted(set(coords[chrom])):
        rows = coords[chrom] == key
        pos = coords["pos"][rows]
        midpoint = float(min(pos[rows]) + max(pos[rows]))/2
        ticks.append(midpoint)
    
    return ticks
