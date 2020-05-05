"""
Utility functions for plotting.
"""

__author__ = "Fabian Neumann (KIT)"
__copyright__ = "Copyright 2019-2020 Fabian Neumann (KIT), GNU GPL 3"

import yaml
from six import iterkeys, itervalues
import pandas as pd


def load_config():
    with open("../config.yaml", "r") as stream:
        return yaml.safe_load(stream)["plotting"]


def get_curtailment(net):
    generated = net.generators_t.p.filter(regex="solar|ror|wind", axis=1).sum().sum()
    possible = (
        (
            net.generators_t.p_max_pu
            * net.generators.filter(regex="solar|ror|wind", axis=0).p_nom_opt
        )
        .sum()
        .sum()
    )
    return (possible - generated) / possible * 100


def get_system_transmission_volume(n, attr=""):
    return (
        (n.lines["s_nom" + attr] * n.lines.length).sum()
        + (n.links["p_nom" + attr] * n.links.length).sum()
    ) / 1e6  # TWkm


def aggregate_costs(n, existing_only=False):

    m = n.copy()

    components = dict(
        Link=("p_nom", "p0"),
        Generator=("p_nom", "p"),
        StorageUnit=("p_nom", "p"),
        Store=("e_nom", "p"),
        Line=("s_nom", None),
    )

    costs = {}
    for c, (p_nom, p_attr) in zip(
        m.iterate_components(iterkeys(components), skip_empty=False),
        itervalues(components),
    ):
        if c.df.empty:
            continue
        if not existing_only:
            p_nom += "_opt"
        if "carrier" not in c.df.columns:
            if c.name == "Line":
                c.df["carrier"] = "AC"
            else:
                c.df["carrier"] = "unspecified"
        costs[(c.list_name, "capital")] = (
            (c.df[p_nom] * c.df.capital_cost).groupby(c.df.carrier).sum()
        )
        if p_attr is not None:
            p = c.pnl[p_attr].sum()
            if c.name == "StorageUnit":
                p = p.loc[p > 0]
            costs[(c.list_name, "marginal")] = (
                (p * c.df.marginal_cost).groupby(c.df.carrier).sum()
            )
    costs = pd.concat(costs)

    return costs
