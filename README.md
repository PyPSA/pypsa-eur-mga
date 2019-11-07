# PyPSA-Eur-MGA: Analysing the Near-Optimal of a Renewable Power System Model

## Paper

Models for long-term investment planning of the
power system typically return a single optimal solution per
set of cost assumptions. However, typically there are many
near-optimal alternatives that stand out due to other attractive
properties like social acceptance. Understanding features that
persist across many cost-efficient alternatives enhances policy
advice and acknowledges structural model uncertainties. We
apply the modeling-to-generate-alternatives (MGA) methodology
to systematically explore the near-optimal feasible space
of a completely renewable European electricity system model.
While accounting for complex spatio-temporal patterns, we allow
simultaneous capacity expansion of generation, storage and
transmission infrastructure subject to linearized multi-period
optimal power flow. Many similarly costly, but technologically
diverse solutions exist. Already a cost deviation of 0.5% offers
a large range of possible investments. However, either offshore
or onshore wind energy along with some hydrogen storage and
transmission network reinforcement are essential to keep costs
within 10% of the optimum.

- F. Neumann, T. Brown, [The Near-Optimal of a Renewable Power System Model](https://arxiv.org/abs/1910.01891), 2019, Preprint submitted to PSCC2020, [arXiv:1910.01891](https://arxiv.org/abs/1910.01891)

## Installation and Usage

Clone the repository

```bash
../ % git clone https://git.scc.kit.edu/FN/pypsa-eur-mga.git
```

and run

```bash
../pypsa-eur-mga % git submodule update --init
```

to retrieve the PyPSA-Eur submodule.

Install and activate the `conda` environment with

```bash
conda env create -f environment.yaml
conda activate pypsa-eur-mga
```

## Workflow

TODO: workflow chart

- `cluster_time`: Clusters network to `{snapshots}` time periods. If `n.storage_units` are present, they are removed since time series clustering distorts order of snapshots.
- `solve_base`: Solves the network to optimality using the regular cost-minimisation objective, which serves as reference value for the MGA iterations.
- `generate_list_of_alternatives`: Generates a list of alternative objectives defining the power system component, its carrier filter, and the optimisation sense as `<COMPONENT>+<CARRIER>+<SENSE>`. Each experiment corresponds to one MGA iteration. This is a [snakemake checkpoint](https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#data-dependent-conditional-execution).
- `generate_alternative`: Solves the network to optimality with the original cost-minimisation objective as constraint with the cost minimum plus some slack `{epsilon}` as upper bound. The new objective is built from the `{objective}` wildcard (from `generate_list_of_alternatives`).
- `extract_results`: Collects and exports results of all near-optimal solutions into several `.csv` files.
- `extract_curtailment`: Extracts curtailment data into a separate `.csv` file.
- `extract_gini`: Calculates and exports the Gini coefficient of different near-optimal solutions.

## Results

TODO: Data on zenodo.
