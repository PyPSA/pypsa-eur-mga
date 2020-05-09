"""
Correlation matrix plotting functions.
"""

__author__ = "Fabian Neumann (KIT)"
__copyright__ = "Copyright 2019-2020 Fabian Neumann (KIT), GNU GPL 3"

import pandas as pd
import numpy as np
import cufflinks as cf
import seaborn as sns
import scipy.cluster.hierarchy as sch
import matplotlib.pyplot as plt

from .utils import load_config

opts = load_config()

cf.go_offline()


def cluster_correlations(df, layered=False):

    cluster_th = 8

    X = df.corr().values
    d = sch.distance.pdist(X)
    L = sch.linkage(d, method="complete")
    ind = sch.fcluster(L, 0.5 * d.max(), "distance")

    columns = [df.columns.tolist()[i] for i in list(np.argsort(ind))]

    if layered:
        unique, counts = np.unique(ind, return_counts=True)
        counts = dict(zip(unique, counts))

        i = 0
        j = 0
        columns = []
        for cluster_l1 in set(sorted(ind)):
            j += counts[cluster_l1]
            sub = df[df.columns.values[i:j]]
            if counts[cluster_l1] > cluster_th:
                X = sub.corr().values
                d = sch.distance.pdist(X)
                L = sch.linkage(d, method="complete")
                ind = sch.fcluster(L, 0.5 * d.max(), "distance")
                col = [sub.columns.tolist()[i] for i in list((np.argsort(ind)))]
                sub = sub.reindex(col, axis="columns")
            cols = sub.columns.tolist()
            columns.extend(cols)
            i = j

    return df.reindex(columns, axis="columns")


def plot_correlation(
    df,
    cluster=False,
    layered_cluster=False,
    iplot=False,
    triangle=False,
    like=None,
    sort_by_carrier=False,
    figsize=(35, 30),
    title="",
):

    if iplot and type(df.columns) == pd.core.indexes.multi.MultiIndex:
        df.columns = [" ".join(col).strip() for col in df.columns.values]

    if cluster:
        df = cluster_correlations(df, layered=layered_cluster)

    if sort_by_carrier:
        columns = []
        for l in [
            "offwind-dc",
            "offwind-ac",
            "onwind",
            "solar",
            "CCGT",
            "OCGT",
            "H2",
            "battery",
            "LK",
            "LN",
        ]:
            columns.extend(df.filter(like=l, axis=1).columns)
        df = df.reindex(columns, axis="columns")

    df = df.rename(columns=opts["nice_names"])
    corr = df.corr()

    if like is not None:
        corr = corr.filter(like=like)

    if iplot:

        if triangle:
            corr = corr.where(np.tril(np.ones(corr.shape)).astype(np.bool))

        lt = cf.Layout(height=1000, width=1000)

        corr.iplot(kind="heatmap", colorscale="PiYG", zmin=-1, zmax=1, layout=lt)

    else:

        f, ax = plt.subplots(figsize=figsize)
        mask = np.zeros_like(corr, dtype=np.bool)

        if triangle:
            mask[np.triu_indices_from(mask)] = True

        sns.heatmap(
            corr,
            mask=mask,
            vmin=-1.0,
            vmax=1.0,
            cmap=sns.diverging_palette(230, 10, as_cmap=True),
            square=True,
            ax=ax,
        )

        ax.set_ylabel("")
        ax.set_xlabel("")
        ax.set_title(title)


def plot_capacity_correlation(
    generation_capacity, storage_capacity, line_volume, link_volume, fn=None
):

    tvol = line_volume.sum() / 1e3 + link_volume.sum() / 1e3
    tvol.name = "Transmission"
    tvol = pd.DataFrame(tvol).T

    selection = [
        False if i in ["PHS", "hydro"] else True
        for i in storage_capacity.index.get_level_values("carrier")
    ]
    sto = storage_capacity.loc[selection].groupby("carrier").sum().sort_index(level=1)

    selection = generation_capacity.index.get_level_values("carrier") != "ror"
    gen = (
        generation_capacity.loc[selection].groupby("carrier").sum().sort_index(level=1)
    )

    df = pd.concat([gen, sto, tvol], axis=0).T

    plot_correlation(
        df,
        triangle=False,
        cluster=True,
        layered_cluster=True,
        figsize=(3.5, 3.5),
        title="Capacity",
    )

    if fn is not None:
        plt.savefig(fn, bbox_inches="tight")


def plot_energy_correlation(energy, fn=None):

    selection = energy.index.get_level_values("carrier") != "ror"
    df = energy.loc[selection].groupby("carrier").sum().sort_index(level=1).T

    plot_correlation(df, triangle=False, figsize=(3.5, 3.5), title="Energy")

    if fn is not None:
        plt.savefig(fn, bbox_inches="tight")
