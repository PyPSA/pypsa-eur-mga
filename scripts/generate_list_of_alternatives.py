import pypsa

n = pypsa.Network(snakemake.input[0])

extendable_components = {
    'generators': n.generators.loc[n.generators.p_nom_extendable].index,
    'storage_units': n.storage_units.loc[n.storage_units.p_nom_extendable].index,
    'stores': n.stores.loc[n.stores.e_nom_extendable].index,
    'links': n.links.loc[n.links.p_nom_extendable].index,
    'lines': n.lines.loc[n.lines.s_nom_extendable].index
}

def mga_components(component_type, component_names):
    return ['+'.join([component_type, c, sense]) for c in component_names for sense in ['min', 'max']]

mga_list = []
for c_type, c_name in extendable_components.items():
    mga_list += mga_components(c_type, c_name)

with open(snakemake.output[0], "w") as f:
    for entry in mga_list:
        f.write(entry +"\n")