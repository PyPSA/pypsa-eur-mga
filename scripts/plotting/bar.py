"""
System cost and capacity bar plotting functions.
"""

__author__ = "Fabian Neumann (KIT)"
__copyright__ = "Copyright 2019-2020 Fabian Neumann (KIT), GNU GPL 3"

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from .utils import load_config, aggregate_costs

opts = load_config()


def plot_cost_bar(n, order=None, relative=True, fn=None):

    tech_colors = opts["tech_colors"]
    nice_names = opts["nice_names"]

    fig, ax = plt.subplots(figsize=(1.5, 5))

    bottom = np.array([0.0, 0.0])

    total_load = (n.snapshot_weightings * n.generators_t.p.sum(axis=1)).sum()

    costs = aggregate_costs(n)
    costs_graph = costs.groupby(level=2).sum()

    if order is not None:
        order = [
            "AC",
            "DC",
            "H2",
            "battery",
            "PHS",
            "hydro",
            "ror",
            "solar",
            "onwind",
            "offwind-ac",
            "offwind-dc",
        ]
        costs_graph = costs_graph.reindex(order)

    for i, ind in enumerate(costs_graph.index):
        data = np.asarray(costs_graph.loc[ind])
        if relative:
            data = data / total_load
        else:
            data = data / 1e9
        ax.bar(
            [0.5],
            data,
            bottom=bottom,
            color=tech_colors[ind],
            linewidth=0,
            width=0.5,
            zorder=-1,
            label=nice_names[ind],
        )
        bottom = bottom + data

    if relative:
        ax.set_ylabel("Average system cost [EUR/MWh]")
        ax.set_ylim([0, 85])
    else:
        ax.set_ylabel("Total system cost [Billion EUR]")
        ax.set_ylim([0, 240])
        ax.set_yticks(np.arange(0, 245, 20))

    ax.set_xlim([0, 1])
    ax.grid(True, axis="y", color="k", linestyle="dotted")
    ax.set_xticks([], [])

    handles, labels = ax.get_legend_handles_labels()
    legend = ax.legend(
        handles[::-1],
        labels[::-1],
        loc="upper left",
        bbox_to_anchor=(1, 1),
        frameon=False,
    )
    ax.add_artist(legend)

    if fn is not None:
        plt.savefig(fn, bbox_inches="tight")


def plot_capacity_bar(n, order=False, fn=None):

    tech_colors = opts["tech_colors"]
    nice_names = opts["nice_names"]

    fig, ax = plt.subplots(figsize=(1.5, 5))

    bottom = np.array([0.0, 0.0])

    gen = n.generators.groupby(n.generators.carrier).p_nom_opt.sum() / 1e3
    sto = n.storage_units.groupby(n.storage_units.carrier).p_nom_opt.sum() / 1e3

    p_nom = pd.concat([gen, sto])

    if order:
        order = [
            "H2",
            "battery",
            "PHS",
            "hydro",
            "ror",
            "solar",
            "onwind",
            "offwind-ac",
            "offwind-dc",
        ]
        p_nom = p_nom.reindex(order)

    for i, ind in enumerate(p_nom.index):
        data = np.asarray(p_nom.loc[ind])
        ax.bar(
            [0.5],
            data,
            bottom=bottom,
            color=tech_colors[ind],
            linewidth=0,
            width=0.5,
            zorder=-1,
            label=nice_names[ind],
        )
        bottom = bottom + data

    ax.set_ylabel("Capacities [GW]")
    ax.set_ylim([0, 1500])
    ax.set_xlim([0, 1])
    ax.grid(True, axis="y", color="k", linestyle="dotted")
    ax.set_xticks([], [])

    handles, labels = ax.get_legend_handles_labels()
    legend = ax.legend(
        handles[::-1],
        labels[::-1],
        loc="upper left",
        bbox_to_anchor=(1, 1),
        frameon=False,
    )
    ax.add_artist(legend)

    if fn is not None:
        plt.savefig(fn, bbox_inches="tight")
