import pypsa
import pandas as pd

def get_mga_components(n):
    extendable_components = {
        'generators': n.generators.loc[n.generators.p_nom_extendable].index,
        'storage_units': n.storage_units.loc[n.storage_units.p_nom_extendable].index,
        'stores': n.stores.loc[n.stores.e_nom_extendable].index,
        'links': n.links.loc[n.links.p_nom_extendable].index,
        'lines': n.lines.loc[n.lines.s_nom_extendable].index
    }
    return extendable_components


def get_mga_groups(n):
    mga_groups = {}
    countries = pd.Series([b[:2] for b in n.buses.index]).unique()
    for comp in ['generators', 'storage_units']:
        carriers = getattr(n, comp).loc[getattr(n, comp).p_nom_extendable].carrier.unique()
        mga_list = list(carriers) + list(countries)
        
        for country in countries:
            country_carriers = getattr(n, comp).loc[
                getattr(n, comp).p_nom_extendable &
                [country in x for x in getattr(n, comp).index]
            ].carrier.unique()
            for carrier in country_carriers:
                mga_list.append(' '.join([country, carrier]))
        mga_groups[comp] = mga_list
    mga_groups['lines'] = [''] # all
    mga_groups['links'] = [''] # all
    mga_groups['transmission'] = [''] # all transmission

    return mga_groups

def mga_list_from_class(component_type, component_names):
    return ['+'.join([component_type, c, sense]) for c in component_names for sense in ['min', 'max']]

n = pypsa.Network(snakemake.input[0])


if snakemake.wildcards.category == "hypercube":
    mga_class = get_mga_components(n) 
elif snakemake.wildcards.category == "groups":
    mga_class = get_mga_groups(n) 
else:
    mga_class = {}

mga_list = []
for c_type, c_name in mga_class.items():
    mga_list += mga_list_from_class(c_type, c_name)

with open(snakemake.output[0], "w") as f:
    for entry in mga_list:
        f.write(entry +"\n")