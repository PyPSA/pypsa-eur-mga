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

from pyomo.environ import Constraint, Objective
from pyomo.core.expr.current import clone_expression

# Add pypsa-eur scripts to path for import
base_path = "/".join(os.getcwd().split("/")[:-1])
sys.path.insert(0, base_path+"/pypsa-eur/scripts")

from solve_network import solve_network, prepare_network, patch_pyomo_tmpdir


def modify_model(n, snapshots, var_type, var_name, obj_sense):


    def encode_objective_as_constraint(n):
        epsilon = float(snakemake.wildcards.epsilon)

        # adding costs of pre-existing infrastructure for suitable base objective value
        constant =  (n.links.p_nom * n.links.capital_cost).sum() + \
            (n.generators.p_nom * n.generators.capital_cost).sum() + \
            (n.lines.s_nom * n.lines.capital_cost).sum() + \
            (n.stores.e_nom * n.stores.capital_cost).sum() + \
            (n.storage_units.p_nom * n.storage_units.capital_cost).sum()
        
        expr = n.model.objective.expr + constant <= (1 + epsilon) * (n.objective + constant)
        n.model.objective_value_slack = Constraint(expr=expr)


    def set_alternative_objective(n, var_types, var_name, obj_sense):

        if type(var_types) is str:
            var_types = [var_types]

        n.model.del_component("objective")

        expr = None
        for var_type in var_types:
            variables = []
            for var in getattr(n.model,var_type):
                
                # line variables are saved with tuple index ('Line', 'var')
                if type(var) is tuple:
                    var_check = var[1]    
                else:
                    var_check = var

                if any(token == var_check for token in var_name.split(' ')):
                    variables.append(var)
                elif all(token in var_check for token in var_name.split(' ')):
                    variables.append(var)

            # translate between pyomo model names and network data names
            m_to_n = {'link_p_nom': 'links',
                    'passive_branch_s_nom': 'lines'}

            if var_type in ['link_p_nom', 'passive_branch_s_nom']:
                index = var[1] if var_type=='passive_branch_s_nom' else var
                part_expr = sum(getattr(n,m_to_n[var_type]).length[index]*getattr(n.model,var_type)[var] for var in variables)
            else:
                part_expr = sum(getattr(n.model,var_type)[var] for var in variables)

            expr = expr + part_expr if expr is not None else part_expr

        n.model.objective = Objective(expr=expr, sense=obj_sense)

        # print objective to console
        n.model.objective.pprint()

    # warmstart
    def set_initial_values(n):

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

    def translate_mga_opts(n, mga_opts):
        network_type_names = ['generators', 'storage_units', 'stores', 'lines', 'links', 'transmission']
        model_type_names = ['generator_p_nom', 'storage_p_nom', 'store_e_nom', 'passive_branch_s_nom', 'link_p_nom', ['passive_branch_s_nom', 'link_p_nom']]
        type_names_dict = dict(zip(network_type_names, model_type_names))
        sense_dict = {'max': -1, 'min': 1}

        def country_pair_component_names(n, country_ids, components):
            if type(components) == str:
                components = [components]
                
            index = []
            for component in components:
                comp_df = getattr(n, component)
                cp = country_ids.split(' ')
                if len(cp) == 1:
                    selector = [(cp[0] == b0[:2] and cp[0] == b1[:2])
                                for b0, b1 in zip(comp_df.bus0, comp_df.bus1)]
                else:
                    selector = [(cp[0] == b0[:2] and cp[1] == b1[:2]) or 
                                (cp[1] == b0[:2] and cp[0] == b1[:2])
                                for b0, b1 in zip(comp_df.bus0, comp_df.bus1)]
                index += list(comp_df.loc[selector].index)

            return " ".join(index)

        subtypes_lookup = {'lines': 'lines', 'links': 'links', 'transmission': ['links', 'lines']}
        if mga_opts[0] in subtypes_lookup.keys():
            mga_opts[1] = country_pair_component_names(n, mga_opts[1], subtypes_lookup[mga_opts[0]])
        mga_opts[0] = type_names_dict[mga_opts[0]]
        mga_opts[2] = sense_dict[mga_opts[2]]

        return mga_opts


    with memory_logger(filename=getattr(snakemake.log, 'memory', None), interval=30.) as mem:
        
        n = pypsa.Network(snakemake.input[0])
        
        mga_opts = translate_mga_opts(n, snakemake.wildcards.objective.split('+'))

        n = prepare_network(n, solve_opts=snakemake.config['solving']['options'])
        
        # catch and tag numerical issues
        try:
            n = solve_network(n, 
                    config=snakemake.config['solving'],
                    solver_log=snakemake.log.solver,
                    opts=opts,
                    extra_functionality=modify_model,
                    extra_functionality_args=mga_opts,
                    skip_iterating=True)
            n.numerical_issue = 0
        except:
            n.numerical_issue = 1

        n.export_to_netcdf(snakemake.output[0])

    logger.info("Maximum memory usage: {}".format(mem.mem_usage))
