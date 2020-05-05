"""
Violin plotting functions.
"""

__author__ = "Fabian Neumann (KIT)"
__copyright__ = "Copyright 2019-2020 Fabian Neumann (KIT), GNU GPL 3"

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

from .utils import load_config
opts = load_config()


def plot_violins(
    generation_capacity, storage_capacity, epsilons=["0.01", "0.05", "0.1"], fn=None
):

    selection = generation_capacity.index.get_level_values("carrier") != "ror"
    gen = generation_capacity.loc[selection] / 1e3

    selection = [
        False if i in ["PHS", "hydro"] else True
        for i in storage_capacity.index.get_level_values("carrier")
    ]
    sto = storage_capacity.loc[selection] / 1e3

    df = pd.concat([gen, sto]).groupby("carrier").sum().unstack().reset_index()

    selection = df.epsilon.apply(lambda e: e in epsilons)
    df = df.loc[selection]

    df.epsilon = df.epsilon.apply(lambda e: float(e) * 100)

    df.carrier = df.carrier.map(opts["nice_names"])

    fig, ax = plt.subplots(figsize=(16, 3))

    sns.violinplot(x=df["carrier"], y=df[0], hue=df["epsilon"], ax=ax)

    ax.set_ylabel("GW")
    ax.set_xlabel("")
    ax.set_ylim([0, 1700])

    plt.legend(title="$\epsilon$ [%]", loc="upper left", fontsize=12)

    if fn is not None:
        plt.savefig(fn, bbox_inches="tight")
