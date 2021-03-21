
import collections
import logging
import numpy as np
import pandas as pd
import toyplot
import toytree

from collections import Counter
from sklearn.linear_model import LinearRegression

from .util import *

def plot_simulation(dat_df,
                    tree,
                    width=700,
                    height=1000,
                    tip_colors=None,
                    log_colors=False,
                    verbose=False):
    canvas = toyplot.Canvas(width=width, height=height)

    ax0 = canvas.cartesian(bounds=('5%', '45%', '5%', '24%'))
    ax1 = canvas.cartesian(bounds=('5%', '45%', '30%', '50%'))
    ax2 = canvas.cartesian(bounds=('5%', '45%', '56%', '74%'))
    ax3 = canvas.cartesian(bounds=('5%', '45%', '80%', '97%'))
    ax4 = canvas.cartesian(bounds=('55%', '90%', '10%', '90%'))

    pis = np.array([x["pi"] for x in dat_df])
    lambdas = np.array([x["lambda_"] for x in dat_df])
    abunds = np.log10(np.array([x["abundance"] for x in dat_df]))
    rs = np.array([x["r"] for x in dat_df])

    _scatter_slope(ax0, pis, lambdas)
    ax0.label.text = "pi/Spec. rate correlation"
    if verbose: print("pi/speciation rate correlation: ", linear_regressor.coef_[0][0])

    _scatter_slope(ax1, rs, lambdas)
    ax1.label.text = "Growth rate/Spec. rate correlation"
    if verbose: print("r/speciation rate correlation: ", linear_regressor.coef_[0][0])

    _scatter_slope(ax2, pis, abunds)
    ax2.label.text = "pi/abundance correlation"
    if verbose: print("pi/abundance correlation: ", linear_regressor.coef_[0][0])

    # add legend
    ax3.plot(sorted(pis/sum(pis), reverse=True))
    ax3.plot(sorted(abunds/sum(abunds), reverse=True))
    ax3.plot(sorted(lambdas/sum(lambdas), reverse=True))
    ax3.label.text = "pi/abund/spec rate distributions"

    tre = toytree.tree(tree)
    tips = tre.get_tip_labels()
    data = None
    try:
        if tip_colors:
            data = np.array([dat_df[t][tip_colors] for t in tips])
    except KeyError:
        raise iBioGenError(BAD_TIP_COLOR_DESCRIPTOR.format(list(dat_df.iloc[0]),
                                                            tip_colors))

    colors = _get_colors(data, log_colors)
    tre.draw(axes=ax4, tip_labels_colors=colors)
    if tip_colors: ax4.label.text = "Tree tips colored by {}".format(tip_colors)

    return canvas


def plot_simulations_summary(params_df,
                             width=700,
                             height=1000,
                             bins=40,
                             plt_zero=False,
                             scale_hills=False,
                             verbose=False):
    canvas = toyplot.Canvas(width=width, height=height)

    pos = params_df[params_df["slope_sign"] == "positive"]
    neg = params_df[params_df["slope_sign"] == "negative"]
    zer = params_df[params_df["slope_sign"] == "zero"]

    # Left column
    ax0 = canvas.cartesian(bounds=('5%', '45%', '3%', '24%'))
    ax1 = canvas.cartesian(bounds=('5%', '45%', '30%', '50%'))
    ax2 = canvas.cartesian(bounds=('5%', '45%', '56%', '74%'))
    ax3 = canvas.cartesian(bounds=('5%', '45%', '80%', '97%'))
    # Right column
    ax4 = canvas.cartesian(bounds=('55%', '95%', '3%', '24%'))
    ax5 = canvas.cartesian(bounds=('55%', '95%', '30%', '50%'))
    ax6 = canvas.cartesian(bounds=('55%', '95%', '56%', '74%'))
    ax7 = canvas.cartesian(bounds=('55%', '95%', '80%', '97%'))

    # Counters can be unordered, but not dict
    cts = Counter(params_df["slope_sign"])
    labs = sorted(cts)
    cts = {lab:cts[lab] for lab in labs}
    if verbose: print("Counts: ", cts)
    ax0.bars(list(cts.values()), title=labs, color=["blue", "red", "black"], opacity=0.5)
    ax0.x.ticks.locator = toyplot.locator.Explicit(labels=labs)
    ax0.label.text = "Counts per significant slope"

    ax1.bars(np.histogram(pos["ClaDS_sigma"], bins), color="red", opacity=0.5)
    ax1.bars(np.histogram(neg["ClaDS_sigma"], bins), color="blue", opacity=0.5)
    ax1.label.text = "ClaDS sigma"

    ax2.bars(np.histogram(pos["turnover_rate"], bins), color="red", opacity=0.5)
    ax2.bars(np.histogram(neg["turnover_rate"], bins), color="blue", opacity=0.5)
    ax2.label.text = "Turnover rate"

    ax4.bars(np.histogram(pos["ClaDS_alpha"], bins), color="red", opacity=0.5)
    ax4.bars(np.histogram(neg["ClaDS_alpha"], bins), color="blue", opacity=0.5)
    ax4.label.text = "ClaDS alpha"

    ax5.scatterplot(pos["ClaDS_alpha"], pos["ClaDS_sigma"], color="red", opacity=0.5)
    ax5.scatterplot(neg["ClaDS_alpha"], neg["ClaDS_sigma"], color="blue", opacity=0.5)
    if plt_zero: ax5.scatterplot(zer["ClaDS_alpha"], zer["ClaDS_sigma"], color="black", size=3, opacity=0.5)
    ax5.label.text = "ClaDS alpha/sigma"

    if scale_hills:
        scaler = params_df["ntaxa"][0]
    else:
        scaler = 1

    ax6.bars(np.histogram(pos["abund_h1"]/scaler, bins), color="red", opacity=0.5)
    ax6.bars(np.histogram(neg["abund_h1"]/scaler, bins), color="blue", opacity=0.5)
    ax6.label.text = "Abundance Hill 1"

    ax3.bars(np.histogram(pos["pi_h1"]/scaler, bins), color="red", opacity=0.5)
    ax3.bars(np.histogram(neg["pi_h1"]/scaler, bins), color="blue", opacity=0.5)
    ax3.label.text = "pi Hill 1"

    ax7.bars(np.histogram(pos["sp_h1"]/scaler, bins), color="red", opacity=0.5)
    ax7.bars(np.histogram(neg["sp_h1"]/scaler, bins), color="blue", opacity=0.5)
    ax7.label.text = "Spec. rate Hill 1"

    return canvas


def plot_random_simulations(params_df,
                            dat_df,
                            width=700,
                            height=1000,
                            upper_lambda=None,
                            marker_colors=None):
    canvas = toyplot.Canvas(width=width, height=height)

    if upper_lambda:
        maxl = np.max(np.array([list(map(lambda y: y["lambda_"], dat_df[x])) for x in dat_df]), axis=0)
        params_df = params_df.iloc[maxl < upper_lambda]
        #dat_df = dat_df.iloc[maxl < upper_lambda]

    pos = params_df[params_df["slope_sign"] == "positive"]
    neg = params_df[params_df["slope_sign"] == "negative"]
    zer = params_df[params_df["slope_sign"] == "zero"]
    
    # Left column
    ax0 = canvas.cartesian(bounds=('5%', '45%', '3%', '24%'))
    ax1 = canvas.cartesian(bounds=('5%', '45%', '30%', '50%'))
    ax2 = canvas.cartesian(bounds=('5%', '45%', '56%', '74%'))
    ax3 = canvas.cartesian(bounds=('5%', '45%', '80%', '97%'))
    # Right column
    ax4 = canvas.cartesian(bounds=('55%', '95%', '3%', '24%'))
    ax5 = canvas.cartesian(bounds=('55%', '95%', '30%', '50%'))
    ax6 = canvas.cartesian(bounds=('55%', '95%', '56%', '74%'))
    ax7 = canvas.cartesian(bounds=('55%', '95%', '80%', '97%'))

    for ax, idx in zip([ax0, ax1, ax2, ax3], np.random.choice(pos.index, 4)):
        pis = np.array([x["pi"] for x in dat_df.iloc[idx]])
        lambdas = np.array([x["lambda_"] for x in dat_df.iloc[idx]])
        abunds = np.log10(np.array([x["abundance"] for x in dat_df.iloc[idx]]))
        rs = np.array([x["r"] for x in dat_df.iloc[idx]])
        _scatter_slope(ax, pis, lambdas, colors=rs)
    ax0.label.text = "Positive slope"

    for ax, idx in zip([ax4, ax5, ax6, ax7], np.random.choice(neg.index, 4)):
        pis = np.array([x["pi"] for x in dat_df.iloc[idx]])
        lambdas = np.array([x["lambda_"] for x in dat_df.iloc[idx]])
        abunds = np.log10(np.array([x["abundance"] for x in dat_df.iloc[idx]]))
        rs = np.array([x["r"] for x in dat_df.iloc[idx]])
        _scatter_slope(ax, pis, lambdas, colors=rs)
    ax4.label.text = "Negative slope"

    return canvas


def _scatter_slope(ax, dat1, dat2, colors='', log_colors=False):
    colors = _get_colors(colors, log_colors)
    ax.scatterplot(dat1, dat2, color=colors)
    linear_regressor = LinearRegression()
    linear_regressor.fit(dat1.reshape(-1, 1), dat2.reshape(-1, 1))
    Y_pred = linear_regressor.predict(dat1.reshape(-1, 1))
    ax.plot(dat1, Y_pred, color='red')


def _get_colors(data, log_colors=False):
    if len(data):
        if log_colors:
            data = np.log10(data)
        minc = min(data)
        maxc = max(data)
        cm = toyplot.color.diverging.map("BlueRed", domain_min=minc, domain_max=maxc)
        colors = cm.colors(data)
    else:
        colors = None

    return colors


# Error messages
BAD_TIP_COLOR_DESCRIPTOR = """
Coloring tree tips must use one of the available attributes:
  {}
You put: {}
"""


if __name__ == "__main__":
    print("Watdo")
