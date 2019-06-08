configfile: "config.yaml"

wildcard_constraints:
    clusters="[0-9]+m?|all",
    ll="(v|c)([0-9\.]+|opt|all)|all",
    opts="[-+a-zA-Z0-9\.]*",
    epsilon="[0-9\.]*",
    #objective="[-+a-zA-Z0-9]*"

subworkflow pypsaeur:
    workdir: "../pypsa-eur"
    snakefile: "../pypsa-eur/Snakefile"
    configfile: "config_pypsaeur.yaml"


# OPTIMAL SOLUTION


rule solve_base:
    input: pypsaeur("networks/elec_s_{clusters}_l{ll}_{opts}.nc")
    output: "results/networks/elec_s_{clusters}_l{ll}_{opts}.nc",
    log:
        solver="logs/elec_s_{clusters}_l{ll}_{opts}_solver.log",
        python="logs/elec_s_{clusters}_l{ll}_{opts}_python.log",
        memory="logs/elec_s_{clusters}_l{ll}_{opts}_memory.log"
    script: "scripts/solve_base.py"

rule solve_all_bases:
    input:
        expand("results/networks/elec_s_{clusters}_l{ll}_{opts}.nc",
                **config['scenario'])

rule summarize_base:



# MODELLING TO GENERATE ALTERNATIVES

def input_generate_all_alternatives():
    wildcards = {**config['scenario'], **config['alternative']}
    input = []
    for cluster in wildcards['clusters']:
        for ll in wildcards['ll']:
            for opts in wildcards['opts']:
                for epsilon in wildcards['epsilon']:
                    alternatives = "results/alternatives/elec_s_{clusters}_l{ll}_{opts}.txt".format(
                        clusters=cluster,
                        ll=ll,
                        opts=opts)
                    obj_list = []
                    with open(alternatives, "r") as f:
                        for line in f:
                            obj_list.append(line.strip())
                    for obj in obj_list:              
                        input.append(
                            "results/networks/elec_s_{clusters}_l{ll}_{opts}_tol{epsilon}_obj-{objective}.nc".format(
                                clusters=cluster,
                                ll=ll,
                                opts=opts,
                                epsilon=epsilon,
                                objective=obj)
                        )
    return input


rule generate_list_of_alternatives:
    input: "results/networks/elec_s_{clusters}_l{ll}_{opts}.nc"
    output: "results/alternatives/elec_s_{clusters}_l{ll}_{opts}.txt"
    script: "scripts/generate_list_of_alternatives.py"

rule generate_alternative:
    input: "results/networks/elec_s_{clusters}_l{ll}_{opts}.nc"
    output: "results/networks/elec_s_{clusters}_l{ll}_{opts}_tol{epsilon}_obj-{objective}.nc"
    log:
        solver="logs/elec_s_{clusters}_l{ll}_{opts}_tol{epsilon}_obj-{objective}_solver.log",
        python="logs/elec_s_{clusters}_l{ll}_{opts}_tol{epsilon}_obj-{objective}_python.log",
        memory="logs/elec_s_{clusters}_l{ll}_{opts}_tol{epsilon}_obj-{objective}_memory.log"
    script: "scripts/generate_alternative.py"

rule generate_all_alternatives:
    input:
        list=expand("results/alternatives/elec_s_{clusters}_l{ll}_{opts}.txt",
            **config['scenario']),
        network=input_generate_all_alternatives()
        


# EVALUATION
