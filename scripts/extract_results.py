import pypsa
import pandas as pd
import re

import logging
logger = logging.getLogger(__name__)

import progressbar as pgb
pgb.streams.wrap_stderr()

logging.basicConfig(level='ERROR')


def infer_wildcards_from_fn(fn):

    d = {}
    l = fn.split('/')[-1].replace('storage_units', 'storageunits').split('_')[2:]
    
    d['clusters'] = l[0]
    d['ll'] = l[1][1:]
    d['snapshots'] = l[2][1:] if l[2][1:] != '' else 'all'
    d['epsilon'] = l[4][3:]
    d['category'] = l[5][4:]
    d['co2'] = re.findall("[0-9]*\.?[0-9]+$", l[3].split('-')[0])[0]
    d['tres'] = l[3].split('-')[1]
    d['type'],d['variable_name'],d['sense'] = l[6][4:-3].split('+')

    if d['variable_name'] == '':
        d['variable_name'] == 'all'

    return d


def get_investments(n, only_extendable=True):
    
    if only_extendable:
        generators_b = n.generators.p_nom_extendable==True
        lines_b = n.lines.s_nom_extendable==True
        stores_b = n.stores.e_nom_extendable==True
        storage_units_b = n.storage_units.p_nom_extendable==True
        links_b = n.links.p_nom_extendable==True
    else:
        generators_b = [True for _ in n.generators.index]
        lines_b = [True for _ in n.lines.index]
        stores_b = [True for _ in n.stores.index]
        storage_units_b = [True for _ in n.storage_units.index]
        links_b = [True for _ in n.links.index]
    
    return pd.concat([
        n.generators[generators_b].p_nom_opt,
        n.lines[lines_b].s_nom_opt,
        n.links[links_b].p_nom_opt,
        n.stores[stores_b].e_nom_opt, 
        n.storage_units[storage_units_b].p_nom_opt])

def get_capacity_mix(component):
    component['country'] = component.bus.apply(lambda x: x[:2])
    return component.groupby(['country', 'carrier']).p_nom_opt.sum()

def get_energy_mix(n):
    n.generators['country'] = n.generators.bus.apply(lambda x: x[:2])
    n.generators['energy'] = n.generators_t.p.multiply(n.snapshot_weightings, axis=0).sum()
    return n.generators.groupby(['country', 'carrier']).energy.sum()

def country_pair(df):
        c0 = df.bus0[:2]
        c1 = df.bus1[:2]
        if c0 == c1:
            return c0
        elif c0 > c1:
            return "{} {}".format(c1, c0)
        else:
            return "{} {}".format(c0, c1)

def get_transmission(n, components, length=False):
    
    if components == 'lines':
        attr = 's_nom_opt'
    else:
        attr = 'p_nom_opt'
        
    getattr(n, components)['country_pair'] = getattr(n, components).apply(country_pair, axis=1)
    
    if length:
        getattr(n, components)['gwkm_opt'] = getattr(n, components)[attr] * getattr(n, components).length / 1e3
        return getattr(n, components).groupby('country_pair').gwkm_opt.sum()
    else:
        return getattr(n, components).groupby('country_pair')[attr].sum()


def get_energy_balance(n, components):
    
    getattr(n, components)['country_pair'] = getattr(n, components).apply(country_pair, axis=1)
    getattr(n, components)['energy'] = getattr(n, components+'_t').p0.multiply(n.snapshot_weightings, axis=0).sum()
    
    def redirect_energy(df):
        c0 = df.bus0[:2]
        c1 = df.bus1[:2]
        if c0 == c1:
            return abs(df.energy)
        elif c0 > c1:
            return -df.energy
        else:
            return df.energy
        
    getattr(n, components)['energy_redirected'] = getattr(n, components).apply(redirect_energy, axis=1)
    return getattr(n, components).groupby('country_pair').energy_redirected.sum()


if __name__ == "__main__":

    ids = ['clusters', 'll', 'snapshots', 'epsilon', 'category', 'co2', 'tres', 'type', 'variable_name', 'sense']
    columns_dict = {n : infer_wildcards_from_fn(n) for n in snakemake.input}
    columns = pd.MultiIndex.from_tuples([e.values() for e in columns_dict.values()], names=ids)
    investment_results = pd.DataFrame(columns=columns)
    gen_capacity_mix = pd.DataFrame(columns=columns)
    sto_capacity_mix = pd.DataFrame(columns=columns)
    energy_mix = pd.DataFrame(columns=columns)
    line_transmission_capacity = pd.DataFrame(columns=columns)
    line_transmission_volume = pd.DataFrame(columns=columns)
    link_transmission_capacity = pd.DataFrame(columns=columns)
    link_transmission_volume = pd.DataFrame(columns=columns)
    line_energy_balance = pd.DataFrame(columns=columns)
    link_energy_balance = pd.DataFrame(columns=columns)

    widgets = [
        pgb.widgets.Percentage(),
        ' ', pgb.widgets.SimpleProgress(format='(%s)' % pgb.widgets.SimpleProgress.DEFAULT_FORMAT),
        ' ', pgb.widgets.Bar(),
        ' ', pgb.widgets.Timer(),
        ' ', pgb.widgets.ETA()
    ]
    progressbar = pgb.ProgressBar(prefix='Extract results: ', widgets=widgets, max_value=len(snakemake.input))

    for k,v in progressbar(columns_dict.items()):

        ids = v.values()
        n = pypsa.Network(k)
        investment_results[(*ids,)] = get_investments(n)
        gen_capacity_mix[(*ids,)] = get_capacity_mix(n.generators)
        sto_capacity_mix[(*ids,)] = get_capacity_mix(n.storage_units)
        energy_mix[(*ids,)] = get_energy_mix(n)
        line_transmission_capacity[(*ids,)] = get_transmission(n, 'lines')
        line_transmission_volume[(*ids,)] = get_transmission(n, 'lines', length=True)
        link_transmission_capacity[(*ids,)] = get_transmission(n, 'links')
        link_transmission_volume[(*ids,)] = get_transmission(n, 'links', length=True)
        line_energy_balance[(*ids,)] = get_energy_balance(n, 'lines')
        link_energy_balance[(*ids,)] = get_energy_balance(n, 'links')

    investment_results.to_csv(snakemake.output.investments)
    energy_mix.to_csv(snakemake.output.energy)
    sto_capacity_mix.to_csv(snakemake.output.storage_capacity)
    gen_capacity_mix.to_csv(snakemake.output.generation_capacity)
    line_transmission_capacity.to_csv(snakemake.output.line_capacity)
    line_transmission_volume.to_csv(snakemake.output.line_volume)
    link_transmission_capacity.to_csv(snakemake.output.link_capacity)
    link_transmission_volume.to_csv(snakemake.output.link_volume)
    line_energy_balance.to_csv(snakemake.output.line_energy_balance)
    link_energy_balance.to_csv(snakemake.output.link_energy_balance)