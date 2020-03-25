import pypsa
import sys
import os
import re

import logging
logger = logging.getLogger(__name__)

# Suppress logging of the slack bus choices
pypsa.pf.logger.setLevel(logging.WARNING)

from vresutils.benchmark import memory_logger

from pyomo.environ import Constraint, Objective
from pyomo.core.expr.current import clone_expression

# Add pypsa-eur scripts to path for import
sys.path.insert(0, os.getcwd() + "/pypsa-eur/scripts")

from solve_network import solve_network, prepare_network, patch_pyomo_tmpdir


def country_pair_component_names(n, country_ids, components):

    if isinstance(components, str):
        components = [components]

    index = []
    for component in components:
        comp_df = getattr(n, component)
        cp = country_ids.split(" ")
        if len(cp) == 1:
            selector = [
                (cp[0] == b0[:2] and cp[0] == b1[:2])
                for b0, b1 in zip(comp_df.bus0, comp_df.bus1)
            ]
        else:
            selector = [
                (cp[0] == b0[:2] and cp[1] == b1[:2])
                or (cp[1] == b0[:2] and cp[0] == b1[:2])
                for b0, b1 in zip(comp_df.bus0, comp_df.bus1)
            ]
        index += list(comp_df.loc[selector].index)

    return " ".join(index)


def translate_mga_opts(n, mga_opts):

        network_type_names = [
            "generators",
            "storage_units",
            "stores",
            "lines",
            "links",
            "transmission",
        ]
        model_type_names = [
            "generator_p_nom",
            "storage_p_nom",
            "store_e_nom",
            "passive_branch_s_nom",
            "link_p_nom",
            ["passive_branch_s_nom", "link_p_nom"],
        ]
        type_names_dict = dict(zip(network_type_names, model_type_names))
        
        sense_dict = {"max": -1, "min": 1}

        subtypes_lookup = {
            "lines": "lines",
            "links": "links",
            "transmission": ["links", "lines"],
        }

        if mga_opts[0] in subtypes_lookup.keys():
            mga_opts[1] = country_pair_component_names(
                n, mga_opts[1], subtypes_lookup[mga_opts[0]]
            )
        mga_opts[0] = type_names_dict[mga_opts[0]]
        mga_opts[2] = sense_dict[mga_opts[2]]

        # attach to network
        n.mga_opts = mga_opts


def encode_objective_as_constraint(n):

    epsilon = float(snakemake.wildcards.epsilon)

    expr = n.model.objective.expr + n.objective_constant <= (1 + epsilon) * (
        n.objective + n.objective_constant
    )

    n.model.objective_value_slack = Constraint(expr=expr)


def set_alternative_objective(n):

    var_type, var_name, obj_sense = n.mga_opts

    if isinstance(var_types, str):
        var_types = [var_types]

    n.model.del_component("objective")

    expr = None
    for var_type in var_types:
        variables = []
        for var in getattr(n.model, var_type):

            # line variables are saved with tuple index ('Line', 'var')
            if isinstance(var, tuple):
                var_check = var[1]
            else:
                var_check = var

            if any(token == var_check for token in var_name.split(" ")):
                variables.append(var)
            elif all(token in var_check for token in var_name.split(" ")):
                variables.append(var)

        # translate between pyomo model names and network data names
        m_to_n = {"link_p_nom": "links", "passive_branch_s_nom": "lines"}

        if var_type in ["link_p_nom", "passive_branch_s_nom"]:
            index = var[1] if var_type == "passive_branch_s_nom" else var
            part_expr = sum(
                getattr(n, m_to_n[var_type]).length[index]
                * getattr(n.model, var_type)[var]
                for var in variables
            )
        else:
            part_expr = sum(getattr(n.model, var_type)[var] for var in variables)

        expr = expr + part_expr if expr is not None else part_expr

    n.model.objective = Objective(expr=expr, sense=obj_sense)

    # print objective to console
    n.model.objective.pprint()


def modify_model(network, snapshots):
    wcs = snakemake.wildcards.objective.split("+")
    translate_mga_opts(network, wcs)
    encode_objective_as_constraint(network)
    set_alternative_objective(network)


if __name__ == "__main__":

    logging.basicConfig(
        filename=snakemake.log.python, level=snakemake.config["logging_level"]
    )

    opts = [
        o
        for o in snakemake.wildcards.opts.split("-")
        if not re.match(r"^\d+h$", o, re.IGNORECASE)
    ]

    solve_opts = snakemake.config["solving"]

    with memory_logger(
        filename=getattr(snakemake.log, "memory", None), interval=30.0
    ) as mem:

        n = pypsa.Network(snakemake.input[0])

        n = prepare_network(n, solve_opts=solve_opts["options"])

        # catch and tag numerical issues
        try:
            n = solve_network(
                n,
                config=solve_opts,
                solver_log=snakemake.log.solver,
                opts=opts,
                extra_functionality=modify_model,
                skip_iterating=True,
            )
            n.numerical_issue = 0
        except:
            n.numerical_issue = 1

        n.export_to_netcdf(snakemake.output[0])

    logger.info("Maximum memory usage: {}".format(mem.mem_usage))
