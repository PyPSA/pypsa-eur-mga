"""
Capacity bar collection plotting functions.
"""

__author__ = "Fabian Neumann (KIT)"
__copyright__ = "Copyright 2019-2020 Fabian Neumann (KIT), GNU GPL 3"

import pandas as pd
import matplotlib.pyplot as plt

from .utils import load_config

opts = load_config()


def subplot_bar_collection(
    generation_capacity,
    storage_capacity,
    ax,
    sense="",
    epsilon="",
    xticks=True,
    title=True,
):

    tech_colors = opts["tech_colors"]
    nice_names = opts["nice_names"]

    techs = [
        "wind",
        "onwind",
        "offwind",
        "solar",
        "H2",
        "battery",
        "all-Transmission",
        "all-StorageUnit",
    ]
    tech_selection = [
        tech in techs
        for tech in generation_capacity.columns.get_level_values("variable_name")
    ]

    def select_sense_epsilon(df, sense, epsilon):
        return (df.columns.get_level_values("sense") == sense) & (
            df.columns.get_level_values("epsilon") == epsilon
        )

    gen = (
        generation_capacity.groupby("carrier")
        .sum()
        .loc[
            :,
            select_sense_epsilon(generation_capacity, sense, epsilon) & tech_selection,
        ]
    )

    sto = (
        storage_capacity.groupby("carrier")
        .sum()
        .loc[:, select_sense_epsilon(storage_capacity, sense, epsilon) & tech_selection]
    )

    df = pd.concat([gen, sto])

    df = (
        df.droplevel(
            [name for name in df.columns.names if name != "variable_name"], axis=1
        )
        / 1e3
    )[techs].T

    colors = [tech_colors[i] for i in df.columns]
    df.rename(columns=nice_names, index=nice_names, inplace=True)

    df.plot.bar(stacked=True, color=colors, legend=False, width=0.95, ax=ax)

    if not xticks:
        ax.set_xticklabels([])

    ax.grid(axis="x")
    ax.set_ylabel(f"{sense.upper()}\nGW")
    ax.set_xlabel("")

    if title:
        ax.set_title(f"$\epsilon$: {float(epsilon)*100}%", size=11)


def plot_bar_collection(generation_capacity, storage_capacity, fn=None):

    fig, ax = plt.subplots(2, 3, figsize=(5.5, 5), sharex=False, sharey=True)

    args = [generation_capacity, storage_capacity]

    top_kwargs = {"sense": "min", "xticks": False}
    subplot_bar_collection(*args, ax[0, 0], epsilon="0.01", **top_kwargs)
    subplot_bar_collection(*args, ax[0, 1], epsilon="0.05", **top_kwargs)
    subplot_bar_collection(*args, ax[0, 2], epsilon="0.1", **top_kwargs)

    bot_kwargs = {"sense": "max", "title": False}
    subplot_bar_collection(*args, ax[1, 0], epsilon="0.01", **bot_kwargs)
    subplot_bar_collection(*args, ax[1, 1], epsilon="0.05", **bot_kwargs)
    subplot_bar_collection(*args, ax[1, 2], epsilon="0.1", **bot_kwargs)

    legend_ax = ax[0, 2]
    handles, labels = legend_ax.get_legend_handles_labels()

    legend = plt.figlegend(
        handles[::-1],
        labels[::-1],
        loc="upper center",
        bbox_to_anchor=(0, 0.15, 1.14, 1),
        bbox_transform=plt.gcf().transFigure,
        handletextpad=0.0,
        columnspacing=0.5,
        ncol=3,
    )
    fig.add_artist(legend)

    plt.tight_layout()

    if fn is not None:
        plt.savefig(fn, bbox_inches="tight")
