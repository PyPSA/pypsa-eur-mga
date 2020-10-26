import pypsa
import sys
import os
import re

from pypsa.descriptors import nominal_attrs
from pypsa.linopt import get_var, linexpr, write_objective, define_constraints
from pypsa.linopf import lookup, network_lopf, ilopf
from pypsa.pf import get_switchable_as_dense as get_as_dense
from pypsa.descriptors import get_extendable_i, get_non_extendable_i

import pandas as pd

import logging

logger = logging.getLogger(__name__)

# Suppress logging of the slack bus choices
pypsa.pf.logger.setLevel(logging.WARNING)

from vresutils.benchmark import memory_logger

# Add pypsa-eur scripts to path for import
sys.path.insert(0, os.getcwd() + "/pypsa-eur/scripts")

from solve_network import prepare_network, add_battery_constraints


def to_regex(pattern):
    """[summary]
    """
    return "(" + ").*(".join(pattern.split(" ")) + ")"


def transmission_countries_to_index(n, countries, components):
    """[summary]
    """

    index = []
    for c in components:
        cntrs = countries.split(" ")
        bus_pairs = zip(n.df(c).bus0, n.df(c).bus1)

        if len(cntrs) == 1:  # within a country
            selection = [
                (cntrs[0] == b0[:2] and cntrs[0] == b1[:2]) for b0, b1 in bus_pairs
            ]
        else:  # across countries
            selection = [
                (cntrs[0] == b0[:2] and cntrs[1] == b1[:2])
                or (cntrs[1] == b0[:2] and cntrs[0] == b1[:2])
                for b0, b1 in bus_pairs
            ]
        index += list(n.df(c).loc[selection].index)

    if len(components) > 1 and countries == "":
        return "LN|LK"
    elif len(index) == 0:
        return ""
    else:
        return "^" + "$|^".join(index) + "$"  # regex for exact match


def objective_constant(n, ext=True, nonext=True):
    """[summary]
    """

    if not (ext or nonext):
        return 0.0

    constant = 0.0
    for c, attr in nominal_attrs.items():
        i = pd.Index([])
        if ext:
            i = i.append(get_extendable_i(n, c))
        if nonext:
            i = i.append(get_non_extendable_i(n, c))
        constant += n.df(c)[attr][i] @ n.df(c).capital_cost[i]

    return constant


def process_objective_wildcard(n, mga_obj):
    """[summary]

    Parameters
    ----------
    n : pypsa.Network
    mga_obj : list-like
        [var_type, pattern, sense]
    """

    lookup = {
        "Line": ["Line"],
        "Transmission": ["Link", "Line"],
    }
    if mga_obj[0] in lookup.keys():
        mga_obj[0] = lookup[mga_obj[0]]
        mga_obj[1] = transmission_countries_to_index(n, mga_obj[1], mga_obj[0])

    lookup = {"max": -1, "min": 1}
    mga_obj[2] = lookup[mga_obj[2]]

    # attach to network
    n.mga_obj = mga_obj

    # print mga_obj to console
    print(mga_obj)


def define_mga_constraint(n, sns, epsilon=None, with_fix=None):
    """Build constraint defining near-optimal feasible space

    Parameters
    ----------
    n : pypsa.Network
    sns : Series|list-like
        snapshots
    epsilon : float, optional
        Allowed added cost compared to least-cost solution, by default None
    with_fix : bool, optional
        Calculation of allowed cost penalty should include cost of non-extendable components, by default None
    """

    if epsilon is None:
        epsilon = float(snakemake.wildcards.epsilon)

    if with_fix is None:
        with_fix = snakemake.config.get("include_non_extendable", True)

    expr = []

    # operation
    for c, attr in lookup.query("marginal_cost").index:
        cost = (
            get_as_dense(n, c, "marginal_cost", sns)
            .loc[:, lambda ds: (ds != 0).all()]
            .mul(n.snapshot_weightings[sns], axis=0)
        )
        if cost.empty:
            continue
        expr.append(linexpr((cost, get_var(n, c, attr).loc[sns, cost.columns])).stack())

    # investment
    for c, attr in nominal_attrs.items():
        cost = n.df(c)["capital_cost"][get_extendable_i(n, c)]
        if cost.empty:
            continue
        expr.append(linexpr((cost, get_var(n, c, attr)[cost.index])))

    lhs = pd.concat(expr).sum()

    if with_fix:
        ext_const = objective_constant(n, ext=True, nonext=False)
        nonext_const = objective_constant(n, ext=False, nonext=True)
        rhs = (1 + epsilon) * (n.objective + ext_const + nonext_const) - nonext_const
    else:
        ext_const = objective_constant(n)
        rhs = (1 + epsilon) * (n.objective + ext_const)

    define_constraints(n, lhs, "<=", rhs, "GlobalConstraint", "mu_epsilon")


def define_mga_objective(n):

    components, pattern, sense = n.mga_obj

    if isinstance(components, str):
        components = [components]

    terms = []
    for c in components:

        variables = get_var(n, c, nominal_attrs[c]).filter(regex=to_regex(pattern))

        if c in ["Link", "Line"] and pattern in ["", "LN|LK", "LK|LN"]:
            coeffs = sense * n.df(c).loc[variables.index, "length"]
        else:
            coeffs = sense

        terms.append(linexpr((coeffs, variables)))

    joint_terms = pd.concat(terms)

    write_objective(n, joint_terms)

    # print objective to console
    print(joint_terms)


def to_mga_model(n, sns):
    """Calls extra functionality modules.
    """
    wc = snakemake.wildcards.objective.split("+")
    process_objective_wildcard(n, wc)
    define_mga_objective(n)
    define_mga_constraint(n, sns)
    add_battery_constraints(n)


# adapted from pypsa-eur/scripts/solve_network.py
def solve_network(
    n, config, solver_log=None, opts="", extra_functionality=None, **kwargs
):
    solver_options = config["solving"]["solver"].copy()
    solver_name = solver_options.pop("name")
    track_iterations = config["solving"]["options"].get("track_iterations", False)
    min_iterations = config["solving"]["options"].get("min_iterations", 4)
    max_iterations = config["solving"]["options"].get("max_iterations", 6)

    # add to network for extra_functionality
    n.config = config
    n.opts = opts

    if config["solving"]["options"].get("skip_iterations", False):
        network_lopf(
            n,
            solver_name=solver_name,
            solver_options=solver_options,
            extra_functionality=extra_functionality,
            **kwargs
        )
    else:
        ilopf(
            n,
            solver_name=solver_name,
            solver_options=solver_options,
            track_iterations=track_iterations,
            min_iterations=min_iterations,
            max_iterations=max_iterations,
            extra_functionality=extra_functionality,
            **kwargs
        )
    return n


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

        n = prepare_network(n, solve_opts=snakemake.config["solving"]["options"])

        # # catch and tag numerical issues
        # try:
        #     n.numerical_issue = 0
        # except:
        #     n.numerical_issue = 1
        n = solve_network(
            n,
            config=snakemake.config,
            solver_log=snakemake.log.solver,
            opts=opts,
            extra_functionality=to_mga_model,
            skip_objective=True,
        )

        n.export_to_netcdf(snakemake.output[0])

    logger.info("Maximum memory usage: {}".format(mem.mem_usage))
