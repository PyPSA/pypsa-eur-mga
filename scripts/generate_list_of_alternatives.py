import pypsa
import pandas as pd
import numpy as np

def get_mga_components(n):
    extendable_components = {
        'generators': n.generators.loc[n.generators.p_nom_extendable].index,
        'storage_units': n.storage_units.loc[n.storage_units.p_nom_extendable].index,
        'stores': n.stores.loc[n.stores.e_nom_extendable].index,
        'links': n.links.loc[n.links.p_nom_extendable].index,
        'lines': n.lines.loc[n.lines.s_nom_extendable].index
    }
    return extendable_components

def append_supertype_carriers(carriers):
    # wind to include onwind, offwind-ac, offwind-dc
    # offwind to include offwind-ac, offwind-dc
    # cgt to include ccgt, ocgt
    supertype_carriers = ["wind", "offwind", "CGT"]
    for sc in supertype_carriers:
        if any([sc in c for c in carriers]):
            carriers.append(sc)
    return carriers

def remove_excluded_carriers(carriers):
    excluded_carriers = snakemake.config['excluded_carriers']
    carriers = [c for c in carriers if c not in excluded_carriers]
    return carriers

def get_mga_components_totals(n):

    mga_totals = {}

    # storage_units and generators
    for comp in ['generators', 'storage_units']:
        carriers = list(getattr(n, comp).loc[getattr(n, comp).p_nom_extendable].carrier.unique())
        carriers = append_supertype_carriers(carriers)
        carriers = remove_excluded_carriers(carriers)
        mga_totals[comp] = carriers

    if len(n.storage_units) > 0:
        mga_totals['storage_units'] += ['']

    # lines and links
    mga_totals['transmission'] = ['']
    if snakemake.config['lines_and_links_separate']:
        mga_totals['lines'] = ['']
        mga_totals['links'] = ['']

    return mga_totals

def get_mga_component_groups(n):
    
    mga_groups = {}
    for comp in ['generators', 'storage_units']:
        countries = list(np.unique([b[:2] for b in getattr(n, comp).bus]))
        carriers = list(getattr(n, comp).loc[getattr(n, comp).p_nom_extendable].carrier.unique())
        carriers = append_supertype_carriers(carriers)    
        
        mga_list = []
        for country in countries:
            country_carriers = list(getattr(n, comp).loc[
                getattr(n, comp).p_nom_extendable &
                [country in x for x in getattr(n, comp).index]
            ].carrier.unique())
            country_carriers = append_supertype_carriers(country_carriers)
            for carrier in country_carriers:
                mga_list.append(' '.join([country, carrier]))

        mga_groups[comp] = mga_list
    
    def country_pair(l):
        # ordering of country pair is unique
        # by sorting alphabetically
        country_a = l.bus0[:2]
        country_b = l.bus1[:2]
        if country_a == country_b:
            return country_b
        elif country_a > country_b:
            return "{} {}".format(country_b, country_a)
        else:
            return "{} {}".format(country_a, country_b)
    

    mga_groups['lines'] = list(n.lines.apply(country_pair, axis=1).unique())
    mga_groups['links'] = list(n.links.apply(country_pair, axis=1).unique())
    mga_groups['transmission'] = list(np.unique(np.concatenate((mga_groups['lines'], mga_groups['links']))))
    if not snakemake.config['lines_and_links_separate']:
        mga_groups['lines'] = []
        mga_groups['links'] = []
    
    return mga_groups

def mga_list_from_class(component_type, component_names):
    return ['+'.join([component_type, c, sense]) for c in component_names for sense in ['min', 'max']]

if __name__ == "__main__":

    n = pypsa.Network(snakemake.input[0])

    if snakemake.wildcards.category == "hypercube":
        mga_class = get_mga_components(n) 
    elif snakemake.wildcards.category == "groups":
        mga_class = get_mga_component_groups(n) 
    elif snakemake.wildcards.category == "totals":
        mga_class = get_mga_components_totals(n)
    else:
        mga_class = {}

    mga_list = []
    for c_type, c_name in mga_class.items():
        mga_list += mga_list_from_class(c_type, c_name)

    with open(snakemake.output[0], "w") as f:
        for entry in mga_list:
            f.write(entry +"\n")