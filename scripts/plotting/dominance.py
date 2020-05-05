"""
Dominance map plotting functions.
"""

__author__ = "Fabian Neumann (KIT)"
__copyright__ = "Copyright 2019-2020 Fabian Neumann (KIT), GNU GPL 3"

import matplotlib.pyplot as plt
import geopandas as gpd
import numpy as np

from .utils import load_config

opts = load_config()


def subplot_dominance(n, ax, shapes=None, attr="energy", threshold=None, title=""):

    if shapes is None:
        shapes = gpd.read_file(snakemake.input.shapes)

    tech_colors = opts["tech_colors"]

    if attr == "energy":
        n.generators[attr] = n.generators_t.p.multiply(
            n.snapshot_weightings, axis=0
        ).sum(axis=0)

    df = n.generators.groupby(["bus", "carrier"])[attr].max().unstack()

    alpha = df.divide(df.sum(axis=1), axis=0).max(axis=1)
    carrier = df.idxmax(axis=1)

    if threshold is not None:
        carrier.loc[df.max(axis=1) <= threshold] = "none"

    shapes["alpha"] = alpha.values
    shapes["carrier"] = carrier.values

    for c in shapes.carrier.unique():
        shapes.loc[shapes.carrier == c].plot(color=tech_colors[c], ax=ax)

    ax.axis("off")

    ax.set_title(f"$\epsilon$ = {float(title)*100}%", fontsize=18)


def plot_dominance(
    n, n_tmin1, n_tmin5, n_tmin10, attr="energy", threshold=None, shapes=None, fn=None
):

    if threshold is None:
        # MWh for energy, MW for capacity
        threshold = 10e6 if attr == "energy" else 5e3

    if shapes is None:
        shapes = gpd.read_file(snakemake.input.shapes)

    tech_colors = opts["tech_colors"]
    nice_names = opts["nice_names"]

    fig, ax = plt.subplots(1, 4, figsize=(26, 6))

    kwargs = {"attr": attr, "threshold": threshold, "shapes": shapes}
    subplot_dominance(n, ax[0], title="0.0", **kwargs)
    subplot_dominance(n_tmin1, ax[1], title="0.01", **kwargs)
    subplot_dominance(n_tmin5, ax[2], title="0.05", **kwargs)
    subplot_dominance(n_tmin10, ax[3], title="0.1", **kwargs)

    handles = []
    labels = []
    for t in np.sort(n.generators.carrier.unique()):
        handles.append(
            plt.Line2D(
                [0], [0], color=tech_colors[t], marker="o", markersize=8, linewidth=0
            )
        )
        labels.append(nice_names[t])

    handles.append(
        plt.Line2D([0], [0], color="gray", marker="o", markersize=8, linewidth=0)
    )
    if attr == "energy":
        labels.append(f"Energy Generated < {threshold/1e6} TWh")
    else:
        labels.append(f"Total Capacity < {threshold/1e3} GW")

    legend_ax = ax[3]
    legend = legend_ax.legend(
        handles,
        labels,
        loc="upper center",
        bbox_to_anchor=(0.4, 1),
        handletextpad=0.0,
        columnspacing=0.5,
        ncol=2,
        fontsize=13,
        framealpha=0.5,
    )
    legend_ax.add_artist(legend)

    plt.subplots_adjust(hspace=-0.1, wspace=-0.1)

    if fn is not None:
        plt.savefig(fn, bbox_inches="tight")
