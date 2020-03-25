import pypsa
import sys
import os
import re

import logging
logger = logging.getLogger(__name__)

from vresutils.benchmark import memory_logger

# Suppress logging of the slack bus choices
pypsa.pf.logger.setLevel(logging.WARNING)

# Add pypsa-eur scripts to path for import
sys.path.insert(0, os.getcwd() + "/pypsa-eur/scripts")

from solve_network import solve_network, prepare_network


def adjust_network(n):

    n.lines.s_max_pu = 0.7

    # need unique naming between links and lines
    # when objective includes lines and links
    n.lines.index = ["LN{}".format(i) for i in n.lines.index]
    n.links.index = ["LK{}".format(i) for i in n.links.index]

    n.lines.s_nom_min = n.lines.s_nom
    n.links.p_nom_min = n.links.p_nom

    n.lines.s_nom_max = n.lines.apply(
        lambda x: x.s_nom
        + max(
            x.s_nom,
            float(snakemake.config["line_additional_circuits"])
            * x.s_nom
            / max(x.num_parallel, 1),
        ),
        axis=1,
    )
    n.links.p_nom_max = float(snakemake.config["link_p_nom_max"])


if __name__ == "__main__":

    logging.basicConfig(
        filename=snakemake.log.python, level=snakemake.config["logging_level"]
    )

    opts = [
        o
        for o in snakemake.wildcards.opts.split("-")
        if not re.match(r"^\d+h$", o, re.IGNORECASE)
    ]

    with memory_logger(
        filename=getattr(snakemake.log, "memory", None), interval=30.0
    ) as mem:
        n = pypsa.Network(snakemake.input[0])

        adjust_network(n)

        n = prepare_network(n, solve_opts=snakemake.config["solving"]["options"])
        n = solve_network(
            n,
            config=snakemake.config["solving"],
            solver_log=snakemake.log.solver,
            opts=opts,
            skip_iterating=True,
        )

        n.export_to_netcdf(snakemake.output[0])

    logger.info("Maximum memory usage: {}".format(mem.mem_usage))
