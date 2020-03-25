import pypsa
import pandas as pd

import progressbar as pgb
pgb.streams.wrap_stderr()

import logging
logger = logging.getLogger(__name__)
logging.basicConfig(level="ERROR")

from extract_results import infer_wildcards_from_fn


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
    curtailment = pd.DataFrame(columns=columns, index=["gini"])

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
        curtailment[(*ids,)] = get_curtailment(n)

    curtailment.to_csv(snakemake.output[0])
