"""
Pie chart plotting functions.
"""

__author__ = "Fabian Neumann (KIT)"
__copyright__ = "Copyright 2019-2020 Fabian Neumann (KIT), GNU GPL 3"

import pandas as pd
import matplotlib.pyplot as plt

from .utils import load_config
opts = load_config()


def energy_by_carrier(n):
    return pd.concat(
        [
            n.generators_t.p.multiply(n.snapshot_weightings, axis=0)
            .sum()
            .groupby(n.generators.carrier)
            .sum()
            .loc[lambda s: s > 0],
            n.storage_units_t.inflow.multiply(n.snapshot_weightings, axis=0)
            .sum()
            .groupby(n.storage_units.carrier)
            .sum(),
        ],
        axis=0,
    )


def plot_energy_pie(n, order=True, fn=None):

    e_primary = energy_by_carrier(n)
    fig, ax = plt.subplots()
    fig.set_size_inches(5, 5)

    if order:
        order = ["offwind-dc", "offwind-ac", "onwind", "solar", "ror", "hydro"]
        e_primary = e_primary.reindex(order)
    else:
        order = list(e_primary.index)

    patches, texts, autotexts = ax.pie(
        e_primary,
        startangle=90,
        colors=pd.Series(order).map(opts["tech_colors"]),
        labels=pd.Series(e_primary.index).map(opts["nice_names"]),
        autopct="%.0f%%",
        shadow=False,
    )

    for t1, t2, i in zip(texts, autotexts, e_primary.index):
        if e_primary.at[i] < 0.005 * e_primary.sum():
            t1.remove()
            t2.remove()

    plt.tight_layout()

    if fn is not None:
        plt.savefig(fn)
