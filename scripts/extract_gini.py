import pypsa
import pandas as pd
import re

import logging

logger = logging.getLogger(__name__)

import progressbar as pgb

pgb.streams.wrap_stderr()

logging.basicConfig(level="ERROR")

from extract_results import infer_wildcards_from_fn


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
    dfc = cumulative_share(n)
    return 1 - dfc.load.diff().multiply(dfc.energy.rolling(2).sum(), axis=0).sum()


if __name__ == "__main__":

    ids = [
        "clusters",
        "ll",
        "snapshots",
        "epsilon",
        "category",
        "co2",
        "tres",
        "type",
        "variable_name",
        "sense",
    ]
    columns_dict = {n: infer_wildcards_from_fn(n) for n in snakemake.input}
    columns = pd.MultiIndex.from_tuples(
        [e.values() for e in columns_dict.values()], names=ids
    )
    gini_results = pd.DataFrame(columns=columns, index=["gini"])

    widgets = [
        pgb.widgets.Percentage(),
        " ",
        pgb.widgets.SimpleProgress(
            format="(%s)" % pgb.widgets.SimpleProgress.DEFAULT_FORMAT
        ),
        " ",
        pgb.widgets.Bar(),
        " ",
        pgb.widgets.Timer(),
        " ",
        pgb.widgets.ETA(),
    ]
    progressbar = pgb.ProgressBar(
        prefix="Extract results: ", widgets=widgets, max_value=len(snakemake.input)
    )

    for k, v in progressbar(columns_dict.items()):

        ids = v.values()
        n = pypsa.Network(k)
        gini_results[(*ids,)] = get_gini(n)

    gini_results.to_csv(snakemake.output[0])
