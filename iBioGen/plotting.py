
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

    ax0.scatterplot(pis, lambdas)
    linear_regressor = LinearRegression()
    linear_regressor.fit(pis.reshape(-1, 1), lambdas.reshape(-1, 1))
    Y_pred = linear_regressor.predict(pis.reshape(-1, 1))
    ax0.plot(pis, Y_pred, color='red')
    ax0.label.text = "pi/Spec. rate correlation"
    if verbose: print("pi/speciation rate correlation: ", linear_regressor.coef_[0][0])

    ax1.scatterplot(rs, lambdas)
    linear_regressor = LinearRegression()
    linear_regressor.fit(rs.reshape(-1, 1), lambdas.reshape(-1, 1))
    Y_pred = linear_regressor.predict(rs.reshape(-1, 1))
    ax1.plot(rs, Y_pred, color='red')
    ax1.label.text = "Growth rate/Spec. rate correlation"
    if verbose: print("r/speciation rate correlation: ", linear_regressor.coef_[0][0])

    ax2.scatterplot(pis, abunds)
    linear_regressor = LinearRegression()
    linear_regressor.fit(pis.reshape(-1, 1), abunds.reshape(-1, 1))
    Y_pred = linear_regressor.predict(pis.reshape(-1, 1))
    ax2.plot(pis, Y_pred, color='red')
    ax2.label.text = "pi/abundance correlation"
    if verbose: print("pi/abundance correlation: ", linear_regressor.coef_[0][0])

    # add legend
    ax3.plot(sorted(pis/sum(pis), reverse=True))
    ax3.plot(sorted(abunds/sum(abunds), reverse=True))
    ax3.plot(sorted(lambdas/sum(lambdas), reverse=True))
    ax3.label.text = "pi/abund/spec rate distributions"

    tre = toytree.tree(tree)
    tre.draw(axes=ax4)

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


if __name__ == "__main__":
    print("Watdo")
