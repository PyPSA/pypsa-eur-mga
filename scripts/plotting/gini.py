"""
Distributional equity plotting functions.
"""

__author__ = "Fabian Neumann (KIT)"
__copyright__ = "Copyright 2019-2020 Fabian Neumann (KIT), GNU GPL 3"

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from .utils import load_config, get_system_transmission_volume

opts = load_config()


def cumulative_share(n, by="bus"):

    n.loads["load"] = n.loads_t.p.multiply(n.snapshot_weightings, axis=0).sum()
    n.generators["energy"] = n.generators_t.p.multiply(
        n.snapshot_weightings, axis=0
    ).sum()

    if by == "country":
        n.loads["country"] = n.loads.apply(lambda x: x.bus[:2], axis=1)
        n.generators["country"] = n.generators.apply(lambda x: x.bus[:2], axis=1)

    energy = n.generators.groupby(by).energy.sum()
    load = n.loads.groupby(by).load.sum()

    df = pd.concat([(energy / energy.sum()), (load / load.sum())], axis=1)
    df.sort_values(by="energy", inplace=True)

    return df.cumsum()


def get_gini(n):
    cumshare = cumulative_share(n)
    return (
        1
        - cumshare.load.diff().multiply(cumshare.energy.rolling(2).sum(), axis=0).sum()
    )


def get_kakwani(n, nbase):
    return get_gini(n) - get_gini(nbase)


def plot_gini(gini, n, annot=None, fn=None):

    fig, ax = plt.subplots(figsize=(2.3, 2.8))
    cb_map = ax.scatter(
        gini.gini, gini.tvol, c=gini.epsilon.astype("float") * 100, alpha=0.6
    )

    ax.scatter(
        get_gini(n),
        get_system_transmission_volume(n, attr="_opt"),
        color="k",
        marker="+",
        s=60,
        label="optimum",
    )
    ax.axhline(
        get_system_transmission_volume(n, "_min"),
        linestyle="--",
        color="k",
        label="today's grid",
    )

    ax.set_xlim([0.2, 0.8])
    ax.set_ylabel("TWkm")
    ax.set_xlabel("Gini [-]")
    ax.legend(frameon=False)

    for k, v in annot.items():
        ax.text(*v, k)

    cb_ax = fig.add_axes([1, 0.1, 0.02, 0.8])
    fig.colorbar(cb_map, cax=cb_ax, label="$\epsilon$ [%]")

    plt.tight_layout()

    if fn is not None:
        plt.savefig(fn, bbox_inches="tight")


def plot_lorentz(networks, fn=None):

    fig, ax = plt.subplots(figsize=(4, 4))

    for label, n in networks.items():
        cumshare = cumulative_share(n)
        plt.plot(cumshare.load, cumshare.energy, label=label)
        plt.plot(np.linspace(0, 1, 100), np.linspace(0, 1, 100), c="k")

    plt.xlabel("Cumulative Share of\n Electricity Demand")
    plt.ylabel("Cumulative Share of\n Electricity Generation")
    plt.legend()

    if fn is not None:
        plt.savefig(fn, bbox_inches="tight")
