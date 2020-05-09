"""
Boxplot plotting functions.
"""

__author__ = "Fabian Neumann (KIT)"
__copyright__ = "Copyright 2019-2020 Fabian Neumann (KIT), GNU GPL 3"

import matplotlib.pyplot as plt

from .utils import load_config

opts = load_config()


def plot_plain_boxplot(df, epsilon="0.01", exclude_fixed=True):
    """df can be generation_capacity, storage_capacity,
    lin{e,k}_{capacity, volume, energy_balance}, energy"""

    if exclude_fixed:
        carrier_selection = [
            c not in ["ror", "PHS", "hydro"]
            for c in df.index.get_level_values("carrier")
        ]
    else:
        carrier_selection = [True for _ in df.index]

    ysize = len(df) / 3

    df.loc[
        carrier_selection, df.columns.get_level_values("epsilon") == "0.01"
    ].T.plot.box(figsize=(30, ysize), whis=[0, 1], vert=False, showfliers=False)

    plt.xticks(rotation=90)


def plot_curated_boxplots(capacities, carrier="onwind", epsilon="0.01", fn=None):

    fig, ax = plt.subplots(figsize=(5, 5))

    df = (
        capacities.loc[
            capacities.index.get_level_values("carrier") == carrier,
            capacities.columns.get_level_values("epsilon") == epsilon,
        ].T.droplevel(1, axis=1)
        / 1e3
    )

    df.plot.box(
        ax=ax,
        whis=[0, 1],
        vert=False,
        showfliers=False,
        color=opts["tech_colors"][carrier],
        widths=0.8,
        notch=True,
        patch_artist=True,
        title=opts["nice_names"][carrier],
    )

    ax.set_xlabel("GW")

    plt.tight_layout()

    if fn is not None:
        plt.savefig(fn, bbox_inches="tight")
