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

from pyomo.environ import Constraint, Objective

def modify_model(n, var_type, var_name, obj_sense):

    def encode_objective_as_constraint(n):
        epsilon = float(snakemake.wildcards.epsilon)
        expr = n.model.objective.expr <= (1 + epsilon) * n.objective
        n.model.objective_value_slack = Constraint(expr=expr)

    def set_alternative_objective(n, var_type, var_name, obj_sense):
        n.model.del_component("objective")
        
        variables = []
        for var in getattr(n.model,var_type):
            if all(token in var for token in var_name.split(' ')):
                variables.append(var)

        m_to_n = {'link_p_nom': 'links',
                  'passive_branch_s_nom': 'lines'}
            
        if var_type in ['link_p_nom', 'passive_branch_s_nom']:
            expr = sum(getattr(n,m_to_n[var_type]).length[var]*getattr(n.model,var_type)[var] for var in variables)
        else:
            expr = sum(getattr(n.model,var_type)[var] for var in variables)
            
        sense = obj_sense
        n.model.objective = Objective(expr=expr, sense=sense)

    def set_initial_values(n):
        set_inital_primal_values(n)
        set_inital_dual_values(n)

    def set_inital_primal_values(n):

        for i in n.generators.loc[n.generators.p_nom_extendable].index:
            n.model.generator_p_nom[i].value = n.generators.loc[i].p_nom_opt
            # for t in n.snapshots:
            #    n.model.generator_p[i, t].value = n.generators_t.p.loc[t, i]

        for i in n.storage_units.loc[n.storage_units.p_nom_extendable].index:
            n.model.storage_p_nom[i].value = n.storage_units.loc[i].p_nom_opt
            # for t in n.snapshots:
            #     p = n.storage_units_t.p.loc[t,i]
            #     if p<0:
            #         n.model.storage_p_dispatch[i,t].value = p
            #         n.model.storage_p_store[i,t].value = 0
            #     else:
            #         n.model.storage_p_dispatch[i,t].value = 0
            #         n.model.storage_p_store[i,t].value = p
            #     n.model.state_of_charge[i,t].value = n.storage_units_t.state_of_charge.loc[t,i]
            #     if "hydro" in i:
            #         n.model.storage_p_spill[i,t].value = n.storage_units_t.spill.loc[t,i]

        for i in n.stores.loc[n.stores.e_nom_extendable].index:
            n.model.store_e_nom[i].value = n.stores.loc[i].e_nom_opt
            # for t in n.snapshots:
            #     n.model.store_p[i,t].value = n.stores_t.p.loc[t,i]
            #     n.model.store_e[i,t].value = n.stores_t.e.loc[t,i]

        for i in n.links.loc[n.links.p_nom_extendable].index:
            n.model.link_p_nom[i].value = n.links.loc[i].p_nom_opt
            # for t in n.snapshots:
            #     n.model.link_p[i,t].value = n.links_t.p0.loc[t,i]

        for i in n.lines.loc[n.lines.s_nom_extendable].index:
            n.model.passive_branch_s_nom['Line',i].value = n.lines.loc[i].s_nom_opt
            # for t in n.snapshots:
            #     n.model.passive_branch_p['Line',i,t].value = n.lines_t.p0.loc[t,i]

    def set_inital_dual_values(n):
        pass
    
    
    encode_objective_as_constraint(n)
    set_alternative_objective(n, var_type, var_name, obj_sense)
    set_initial_values(n)

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

    def translate_mga_opts(mga_opts):

        network_type_names = ['generators', 'storage_units', 'stores', 'links', 'lines']
        model_type_names = ['generator_p_nom', 'storage_p_nom', 'store_e_nom', 'passive_branch_s_nom', 'link_p_nom']
        type_names_dict = dict(zip(network_type_names, model_type_names))
        sense_dict = {'max': -1, 'min': 1}

        mga_opts[0] = type_names_dict[mga_opts[0]]
        mga_opts[2] = sense_dict[mga_opts[2]]

        return mga_opts

    mga_opts = translate_mga_opts(snakemake.wildcards.objective.split('+'))

    with memory_logger(filename=getattr(snakemake.log, 'memory', None), interval=30.) as mem:
        n = pypsa.Network(snakemake.input[0])

        n = prepare_network(n, solve_opts=snakemake.config['solving']['options'])
        n = solve_network(n, 
                config=snakemake.config['solving'],
                solver_log=snakemake.log.solver,
                opts=opts,
                modification=modify_model,
                modification_args=mga_opts,
                skip_iterating=True)

        n.export_to_netcdf(snakemake.output[0])

    logger.info("Maximum memory usage: {}".format(mem.mem_usage))
