from itertools import product

configfile: "config.yaml"

wildcard_constraints:
    clusters="[0-9]+m?|all",
    opts="[-+a-zA-Z0-9\.]*",
    epsilon="[0-9\.]*",

subworkflow pypsaeur:
    workdir: "pypsa-eur"
    snakefile: "pypsa-eur/Snakefile"
    configfile: "config.pypsaeur.yaml"


def memory(w):
    factor = 1.3
    for o in w.opts.split('-'):
        m = re.match(r'^(\d+)h$', o, re.IGNORECASE)
        if m is not None:
            factor /= int(m.group(1))
            break
    if w.clusters.endswith('m'):
        return int(factor * (18000 + 180 * int(w.clusters[:-1])))
    else:
        return int(factor * (10000 + 195 * int(w.clusters)))

# OPTIMAL SOLUTION

rule solve_base:
    input: pypsaeur("networks/elec_s_{clusters}_ec_lcopt_{opts}.nc")
    output: "results/networks/elec_s_{clusters}_ec_lcopt_{opts}.nc",
    benchmark: "logs/elec_s_{clusters}_ec_lcopt_{opts}_time.log"
    log:
        solver="logs/elec_s_{clusters}_ec_lcopt_{opts}_solver.log",
        python="logs/elec_s_{clusters}_ec_lcopt_{opts}_python.log",
        memory="logs/elec_s_{clusters}_ec_lcopt_{opts}_memory.log"
    threads: 2
    resources: mem=memory
    script: "scripts/solve_base.py"

rule solve_all_bases:
    input:
        expand("results/networks/elec_s_{clusters}_ec_lcopt_{opts}.nc",
                **config['scenario-totals'])


# MODELLING TO GENERATE ALTERNATIVES

# At this checkpoint (https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#data-dependent-conditional-execution) 
# based on the variables of the original problem the search directions
# of the MGA iterations are inferred.

checkpoint generate_list_of_alternatives:
    input: "results/networks/elec_s_{clusters}_ec_lcopt_{opts}.nc"
    output: "results/alternatives/elec_s_{clusters}_ec_lcopt_{opts}_cat-{category}.txt"
    script: "scripts/generate_list_of_alternatives.py"

rule generate_alternative:
    input: "results/networks/elec_s_{clusters}_ec_lcopt_{opts}.nc"
    output: "results/networks/elec_s_{clusters}_ec_lcopt_{opts}_tol{epsilon}_cat-{category}_obj-{objective}.nc"
    benchmark: "logs/elec_s_{clusters}_ec_lcopt_{opts}_tol{epsilon}_cat-{category}_obj-{objective}_time.log"
    log:
        solver="logs/elec_s_{clusters}_ec_lcopt_{opts}_tol{epsilon}_cat-{category}_obj-{objective}_solver.log",
        python="logs/elec_s_{clusters}_ec_lcopt_{opts}_tol{epsilon}_cat-{category}_obj-{objective}_python.log",
        memory="logs/elec_s_{clusters}_ec_lcopt_{opts}_tol{epsilon}_cat-{category}_obj-{objective}_memory.log"
    threads: 2
    resources: mem=memory
    script: "scripts/generate_alternative.py"


def get_wildcard_sets(config):
    wildcard_sets = [
        {**config['scenario-totals'], **config['alternative-totals']}
    ]
    if config['include_groups']:
        wildcard_sets.append(
            {**config['scenario-groups'], **config['alternative-groups']}
        )
    if config['include_hypercube']:
        wildcard_sets.append(
            {**config['scenario-hypercube'], **config['alternative-hypercube']}
        )
    return wildcard_sets


def append_cost_uncertainty_opts(wildcards_opts):
    factors = config["cost-uncertainty"]
    carrier, values = zip(*factors.items())
    cost_sets = [dict(zip(carrier, v)) for v in product(*values)]

    new_opts = []
    for opts in wildcards_opts:
        for cost_set in cost_sets:
            cost_opts = "-".join([f"{c}+{v}" for c, v in cost_set.items()])
            new_opts.append(f"{opts}-{cost_opts}")
    
    return new_opts


def input_generate_clusters_alternatives(w):
    wildcard_sets = get_wildcard_sets(config)
    input = []
    for wildcards in wildcard_sets:
        for clusters in wildcards["clusters"]:
            if int(clusters) == int(w.clusters):
                wildcards_opts = wildcards["opts"]
                if "cost-uncertainty" in config.keys():
                    wildcards_opts = append_cost_uncertainty_opts(wildcards_opts)
                for opts in wildcards_opts:
                    for epsilon in wildcards['epsilon']:
                        for category in wildcards['category']:
                            alternatives = checkpoints.generate_list_of_alternatives.get(
                                clusters=w.clusters,
                                opts=opts,
                                category=category).output[0]
                            obj_list = []
                            with open(alternatives, "r") as f:  
                                for line in f:
                                    obj_list.append(line.strip())
                            for obj in obj_list:              
                                input.append(
                                    "results/networks/elec_s_{clusters}_ec_lcopt_{opts}_tol{epsilon}_cat-{category}_obj-{objective}.nc".format(
                                        clusters=w.clusters,
                                        opts=opts,
                                        epsilon=epsilon,
                                        objective=obj,
                                        category=category)
                                )
    return input


def input_generate_all_alternatives(w):
    categories = ["totals"]
    if config["include_groups"]: categories.append("groups")
    if config["include_hypercube"]: categories.append("hypercube")
    all_clusters = set().union(*[config[f"scenario-{c}"]["clusters"] for c in categories])
    input = []
    for clusters in all_clusters:
        wcs = snakemake.io.Wildcards(fromdict={"clusters": clusters})
        input.extend(
            input_generate_clusters_alternatives(wcs)
        )
    return input


rule generate_all_alternatives:
    input: input_generate_all_alternatives


# EVALUATION

rule extract_results:
    input: input_generate_clusters_alternatives
    output:
        investments="results/summaries/{clusters}/investments.csv",
        energy="results/summaries/{clusters}/energy.csv",
        storage_capacity="results/summaries/{clusters}/storage_capacity.csv",
        generation_capacity="results/summaries/{clusters}/generation_capacity.csv",
        line_capacity="results/summaries/{clusters}/line_capacity.csv",
        link_capacity="results/summaries/{clusters}/link_capacity.csv",
        line_volume="results/summaries/{clusters}/line_volume.csv",
        link_volume="results/summaries/{clusters}/link_volume.csv",
        line_energy_balance="results/summaries/{clusters}/line_energy_balance.csv",
        link_energy_balance="results/summaries/{clusters}/link_energy_balance.csv",
        curtailment="results/summaries/{clusters}/curtailment.csv",
        gini="results/summaries/{clusters}/gini.csv",
        maps=directory("graphics/{clusters}/networks")
    script: "scripts/extract_results.py"


rule extract_all_results:
    input: expand("results/summaries/{clusters}/investments.csv", clusters=config["scenario-totals"]["clusters"])
