import sys
import os
import numpy as np
import pandas as pd
import logging
logger = logging.getLogger(__name__)
import gc
import re

import pypsa
from pypsa.descriptors import free_output_series_dataframes

# Suppress logging of the slack bus choices
pypsa.pf.logger.setLevel(logging.WARNING)

from vresutils.benchmark import memory_logger

# Add pypsa-eur scripts to path for import
base_path = "/".join(os.getcwd().split("/")[:-1])
sys.path.insert(0, base_path+"/pypsa-eur/scripts")

from solve_network import solve_network, prepare_network, patch_pyomo_tmpdir

if __name__ == "__main__":
    # Detect running outside of snakemake and mock snakemake for testing
    if 'snakemake' not in globals():
        from vresutils.snakemake import MockSnakemake, Dict
        snakemake = MockSnakemake(
            wildcards=dict(network='elec', simpl='', clusters='45', lv='1.0', opts='Co2L-3H'),
            input=["networks/{network}_s{simpl}_{clusters}_lv{lv}_{opts}.nc"],
            output=["results/networks/s{simpl}_{clusters}_lv{lv}_{opts}.nc"],
            log=dict(solver="logs/{network}_s{simpl}_{clusters}_lv{lv}_{opts}_solver.log",
                     python="logs/{network}_s{simpl}_{clusters}_lv{lv}_{opts}_python.log")
        )

    tmpdir = snakemake.config['solving'].get('tmpdir')
    if tmpdir is not None:
        patch_pyomo_tmpdir(tmpdir)

    logging.basicConfig(filename=snakemake.log.python,
                        level=snakemake.config['logging_level'])

    opts = [o
            for o in snakemake.wildcards.opts.split('-')
            if not re.match(r'^\d+h$', o, re.IGNORECASE)]

    with memory_logger(filename=getattr(snakemake.log, 'memory', None), interval=30.) as mem:
        n = pypsa.Network(snakemake.input[0])
        n.lines.s_max_pu = 0.7 # temporary assert
        n.lines.s_nom_min = n.lines.s_nom

        n = prepare_network(n, solve_opts=snakemake.config['solving']['options'])
        n = solve_network(n, config=snakemake.config['solving'], solver_log=snakemake.log.solver, opts=opts, skip_iterating=True)

        n.export_to_netcdf(snakemake.output[0])

    logger.info("Maximum memory usage: {}".format(mem.mem_usage))
