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
    factor = 1.5
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

# generates list of inputs from alternatives list
def input_generate_all_alternatives(w):
    wildcards_sets = [
        {**config['scenario-totals'], **config['alternative-totals']}
    ]
    if config['include_groups']:
        wildcards_sets.append(
            {**config['scenario-groups'], **config['alternative-groups']}
        )
    if config['include_hypercube']:
        wildcards_sets.append(
            {**config['scenario-hypercube'], **config['alternative-hypercube']}
        )
    input = []
    for wildcards in wildcards_sets:
        for cluster in wildcards['clusters']:
            for opts in wildcards['opts']:
                for epsilon in wildcards['epsilon']:
                    for category in wildcards['category']:
                        alternatives = checkpoints.generate_list_of_alternatives.get(
                            clusters=cluster,
                            opts=opts,
                            category=category).output[0]
                        obj_list = []
                        with open(alternatives, "r") as f:  
                            for line in f:
                                obj_list.append(line.strip())
                        for obj in obj_list:              
                            input.append(
                                "results/networks/elec_s_{clusters}_ec_lcopt_{opts}_tol{epsilon}_cat-{category}_obj-{objective}.nc".format(
                                    clusters=cluster,
                                    opts=opts,
                                    epsilon=epsilon,
                                    objective=obj,
                                    category=category)
                            )
    return input

rule generate_all_alternatives:
    input: input_generate_all_alternatives


# EVALUATION

rule extract_results:
    input: input_generate_all_alternatives
    output:
        investments="results/summaries/investments.csv",
        energy="results/summaries/energy.csv",
        storage_capacity="results/summaries/storage_capacity.csv",
        generation_capacity="results/summaries/generation_capacity.csv",
        line_capacity="results/summaries/line_capacity.csv",
        link_capacity="results/summaries/link_capacity.csv",
        line_volume="results/summaries/line_volume.csv",
        link_volume="results/summaries/link_volume.csv",
        line_energy_balance="results/summaries/line_energy_balance.csv",
        link_energy_balance="results/summaries/link_energy_balance.csv"
    script: "scripts/extract_results.py"
    
rule extract_gini:
    input: input_generate_all_alternatives
    output: "results/summaries/gini.csv"
    script: "scripts/extract_gini.py"
    
rule extract_curtailment:
    input: input_generate_all_alternatives
    output: "results/summaries/curtailment.csv"
    script: "scripts/extract_curtailment.py"
    