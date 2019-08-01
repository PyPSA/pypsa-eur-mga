# TODO: fix environment file
import os
os.system("pip install tsam")

configfile: "config.yaml"

wildcard_constraints:
    clusters="[0-9]+m?|all",
    ll="(v|c)([0-9\.]+|opt|all)|all",
    opts="[-+a-zA-Z0-9\.]*",
    epsilon="[0-9\.]*",
    snapshots="[0-9]*"

subworkflow pypsaeur:
    workdir: "../pypsa-eur"
    snakefile: "../pypsa-eur/Snakefile"
    configfile: "config_pypsaeur.yaml"


# OPTIMAL SOLUTION
rule cluster_time:
    input: pypsaeur("networks/elec_s_{clusters}_l{ll}_{opts}.nc")
    output: "networks/elec_s_{clusters}_l{ll}_t{snapshots}_{opts}.nc"
    script: "scripts/cluster_time.py"

rule cluster_all_times:
    input:
        expand("networks/elec_s_{clusters}_l{ll}_t{snapshots}_{opts}.nc",
                **config['scenario'])

rule solve_base:
    input: "networks/elec_s_{clusters}_l{ll}_t{snapshots}_{opts}.nc"
    output: "results/networks/elec_s_{clusters}_l{ll}_t{snapshots}_{opts}.nc",
    benchmark: "logs/elec_s_{clusters}_l{ll}_t{snapshots}_{opts}_time.log"
    log:
        solver="logs/elec_s_{clusters}_l{ll}_t{snapshots}_{opts}_solver.log",
        python="logs/elec_s_{clusters}_l{ll}_t{snapshots}_{opts}_python.log",
        memory="logs/elec_s_{clusters}_l{ll}_t{snapshots}_{opts}_memory.log"
    script: "scripts/solve_base.py"

rule solve_all_bases:
    input:
        expand("results/networks/elec_s_{clusters}_l{ll}_t{snapshots}_{opts}.nc",
                **config['scenario'])


# MODELLING TO GENERATE ALTERNATIVES


checkpoint generate_list_of_alternatives:
    input: "results/networks/elec_s_{clusters}_l{ll}_t{snapshots}_{opts}.nc"
    output: "results/alternatives/elec_s_{clusters}_l{ll}_t{snapshots}_{opts}_cat-{category}.txt"
    script: "scripts/generate_list_of_alternatives.py"

rule generate_alternative:
    input: "results/networks/elec_s_{clusters}_l{ll}_t{snapshots}_{opts}.nc"
    output: "results/networks/elec_s_{clusters}_l{ll}_t{snapshots}_{opts}_tol{epsilon}_cat-{category}_obj-{objective}.nc"
    benchmark: "logs/elec_s_{clusters}_l{ll}_t{snapshots}_{opts}_tol{epsilon}_cat-{category}_obj-{objective}_time.log"
    log:
        solver="logs/elec_s_{clusters}_l{ll}_t{snapshots}_{opts}_tol{epsilon}_cat-{category}_obj-{objective}_solver.log",
        python="logs/elec_s_{clusters}_l{ll}_t{snapshots}_{opts}_tol{epsilon}_cat-{category}_obj-{objective}_python.log",
        memory="logs/elec_s_{clusters}_l{ll}_t{snapshots}_{opts}_tol{epsilon}_cat-{category}_obj-{objective}_memory.log"
    script: "scripts/generate_alternative.py"

def input_generate_all_alternatives(w):
    wildcards = {**config['scenario'], **config['alternative']}
    input = []
    for cluster in wildcards['clusters']:
        for ll in wildcards['ll']:
            for snapshots in wildcards['snapshots']:
                for opts in wildcards['opts']:
                    for epsilon in wildcards['epsilon']:
                        for category in wildcards['category']:
                            alternatives = checkpoints.generate_list_of_alternatives.get(
                                clusters=cluster,
                                ll=ll,
                                snapshots=snapshots,
                                opts=opts,
                                category=category).output[0]
                            obj_list = []
                            with open(alternatives, "r") as f:  
                                for line in f:
                                    obj_list.append(line.strip())
                            for obj in obj_list:              
                                input.append(
                                    "results/networks/elec_s_{clusters}_l{ll}_t{snapshots}_{opts}_tol{epsilon}_cat-{category}_obj-{objective}.nc".format(
                                        clusters=cluster,
                                        ll=ll,
                                        snapshots=snapshots,
                                        opts=opts,
                                        epsilon=epsilon,
                                        objective=obj,
                                        category=category)
                                )
    return input

rule generate_all_alternatives:
    input: input_generate_all_alternatives



# EVALUATION
