"""
Near-optimal system boundaries plotting functions.
"""

__author__ = "Fabian Neumann (KIT)"
__copyright__ = "Copyright 2019-2020 Fabian Neumann (KIT), GNU GPL 3"

import numpy as np
import matplotlib.pyplot as plt

from math import isclose

from .utils import load_config, get_system_transmission_volume

opts = load_config()


def oneport_space(n_opt, investments, c, regex=None):

    if regex is None:
        regex = c

    if regex == "all-StorageUnit":
        regex = "battery|H2"

    if regex in ["battery", "H2", "battery|H2"]:
        comp = "storage_units"
    else:
        comp = "generators"

    selection = investments.columns.get_level_values("variable_name") == c
    to_drop = ["clusters", "variable_name", "category", "type", "tres"]
    df = (
        investments.loc[:, selection]
        .filter(regex=regex, axis=0)
        .sum()
        .droplevel(to_drop)
        / 1e3
    )

    p_nom_opt = getattr(n_opt, comp).filter(regex=regex, axis=0).p_nom_opt.sum() / 1e3

    df.loc[df.apply(lambda i: isclose(i, p_nom_opt))] = np.nan

    df[("0.0", "min")] = p_nom_opt
    df[("0.0", "max")] = p_nom_opt

    if regex == "battery|H2":
        df += (
            n_opt.storage_units.filter(regex="PHS|hydro", axis=0).p_nom_opt.sum() / 1e3
        )

    return df.unstack().reset_index().sort_values(by="epsilon")


def branch_space(n_opt, line_volume, link_volume, c):

    df = link_volume.sum() + line_volume.sum()

    selection = df.index.get_level_values("variable_name") == c
    to_drop = ["clusters", "variable_name", "category", "type", "tres"]
    df = df.loc[:, selection].droplevel(to_drop) / 1e3

    tvol_opt = get_system_transmission_volume(n_opt, "_opt")

    df.loc[df.apply(lambda i: isclose(i, tvol_opt))] = np.nan

    df[("0.0", "min")] = tvol_opt
    df[("0.0", "max")] = tvol_opt

    return df.unstack().reset_index().sort_values(by="epsilon")


def subplot_space(
    n_opt,
    investments,
    line_volume,
    link_volume,
    ax,
    carrier,
    other_carrier=None,
    other_carrier_sense="max",
    xlabel=False,
    ylabel=False,
):

    tech_colors = opts["tech_colors"]
    nice_names = opts["nice_names"]

    if carrier == "all-Transmission":
        df = branch_space(n_opt, line_volume, link_volume, carrier)
        ax.set_ylabel(nice_names[carrier] + " [TWkm]")
    else:
        df = oneport_space(n_opt, investments, carrier)
        ax.set_ylabel(nice_names[carrier] + " [GW]")

    eps = df.epsilon.astype("float") * 100

    kwargs = {"marker": ".", "markersize": 10, "color": tech_colors[carrier]}
    ax.plot(eps, df["max"], **kwargs)
    ax.plot(eps, df["min"], **kwargs)
    ax.fill_between(eps, df["min"], df["max"], alpha=0.5, color=tech_colors[carrier])

    if xlabel:
        ax.set_xlabel("$\epsilon$ [%]")

    if ylabel:
        ax.yaxis.set_tick_params(which="both", labelbottom=True)

    ax.set_xticks(np.arange(0, 11, 2))

    if other_carrier is not None:
        if carrier == "all-Transmission":
            df = branch_space(n_opt, line_volume, link_volume, other_carrier)
        else:
            df = oneport_space(n_opt, investments, other_carrier, carrier)
        eps = df.epsilon.astype("float") * 100
        kwargs = {"color": "k", "linestyle": "-"}
        ax.plot(eps, df[other_carrier_sense], **kwargs)

    if carrier == "all-Transmission":
        ax.axhline(
            get_system_transmission_volume(n_opt, "_min"),
            linestyle="--",
            color="k",
            label="today's grid",
        )
        ax.legend(frameon=False)

    if carrier == "all-StorageUnit":
        ax.axhline(
            n_opt.storage_units.p_nom.sum() / 1e3,
            linestyle="--",
            color="k",
            label="today's\nhydro+PHS",
        )
        ax.legend(frameon=False)


def plot_space(
    n_opt,
    investments,
    line_volume,
    link_volume,
    other_carrier=None,
    other_carrier_sense="max",
    compact=False,
    fn=None,
):

    figsize = (5, 9) if compact else (6, 10)
    fix, ax = plt.subplots(4, 2, figsize=figsize, sharey=True, sharex=compact)

    args = [n_opt, investments, line_volume, link_volume]
    kwargs = {
        "xlabel": not compact,
        "ylabel": not compact,
    }
    oc_kwargs = {
        "other_carrier": other_carrier,
        "other_carrier_sense": other_carrier_sense,
    }
    subplot_space(*args, ax[0, 0], "wind", **kwargs, **oc_kwargs)
    subplot_space(*args, ax[0, 1], "solar", **kwargs, **oc_kwargs)
    subplot_space(*args, ax[1, 0], "onwind", **kwargs, **oc_kwargs)
    subplot_space(*args, ax[1, 1], "offwind", **kwargs, **oc_kwargs)
    subplot_space(*args, ax[2, 0], "all-StorageUnit", **kwargs, **oc_kwargs)
    subplot_space(*args, ax[2, 1], "all-Transmission", **kwargs, **oc_kwargs)
    subplot_space(*args, ax[3, 0], "H2", xlabel=True, ylabel=not compact, **oc_kwargs)
    subplot_space(
        *args, ax[3, 1], "battery", xlabel=True, ylabel=not compact, **oc_kwargs
    )

    plt.tight_layout()

    hspace = 0.1 if compact else 0.35
    plt.subplots_adjust(hspace=hspace)
    
    plt.ylim([-50,1300])

    if fn is not None:
        plt.savefig(fn, bbox_inches="tight")


def plot_space_presentation(
    n_opt,
    investments,
    line_volume,
    link_volume,
    other_carrier=None,
    other_carrier_sense="max",
    fn=None,
):

    figsize = (11.5, 6)
    fix, ax = plt.subplots(2, 4, figsize=figsize, sharey=True, sharex=False)

    args = [n_opt, investments, line_volume, link_volume]
    kwargs = {
        "xlabel": True,
        "ylabel": True,
        "other_carrier": other_carrier,
        "other_carrier_sense": other_carrier_sense,
    }
    subplot_space(*args, ax[0, 0], "wind", **kwargs)
    subplot_space(*args, ax[1, 0], "solar", **kwargs)
    subplot_space(*args, ax[0, 1], "onwind", **kwargs)
    subplot_space(*args, ax[1, 1], "offwind", **kwargs)
    subplot_space(*args, ax[0, 2], "all-StorageUnit", **kwargs)
    subplot_space(*args, ax[1, 2], "all-Transmission", **kwargs)
    subplot_space(*args, ax[0, 3], "H2", **kwargs)
    subplot_space(*args, ax[1, 3], "battery", **kwargs)

    plt.tight_layout()
    
    plt.ylim([-50,1300])

    if fn is not None:
        plt.savefig(fn, bbox_inches="tight")
