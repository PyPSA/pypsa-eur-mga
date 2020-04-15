import pypsa
import numpy as np


def country_pair(branch):
    """Creates a country pair inferred from
    bus names of a line or link.

    Parameters
    ----------
    branch : pd.Series|dict

    Returns
    -------
    str
        Country pair, e.g. "DE FR". Ordering
        is unique since sorted alphabetically.
    """
    country_a = branch.bus0[:2]
    country_b = branch.bus1[:2]
    if country_a == country_b:
        return country_b
    elif country_a > country_b:
        return "{} {}".format(country_b, country_a)
    else:
        return "{} {}".format(country_a, country_b)


def append_supertype_carriers(carriers):
    """Adds supertype carriers to the list of carriers. These
    unite multiple individual carriers if one of its constituents
    is already in the list of carriers. E.g.
    'wind' to include {'onwind', 'offwind-ac', 'offwind-dc'},
    'offwind' to include {'offwind-ac', 'offwind-dc'}, or
    'cgt' to include {'ccgt', 'ocgt'}
    """
    supertype_carriers = ["wind", "offwind", "CGT"]
    for sc in supertype_carriers:
        if any([sc in c for c in carriers]):
            carriers.append(sc)
    return carriers


def remove_excluded_carriers(carriers):
    """Removes carriers listed in the 'config.yaml`
    from the list of carriers.
    """    
    excluded_carriers = snakemake.config["excluded_carriers"]
    return [c for c in carriers if c not in excluded_carriers]


def create_mga_lookup_hypercube(n):
    return {
        "Generator": n.generators.loc[n.generators.p_nom_extendable].index,
        "StorageUnit": n.storage_units.loc[n.storage_units.p_nom_extendable].index,
        "Store": n.stores.loc[n.stores.e_nom_extendable].index,
        "Link": n.links.loc[n.links.p_nom_extendable].index,
        "Line": n.lines.loc[n.lines.s_nom_extendable].index,
    }


def create_mga_lookup_totals(n):

    lookup = {}

    # Storage Units and Generators
    for c in ["Generator", "StorageUnit"]:
        carriers = list(n.df(c).loc[n.df(c).p_nom_extendable].carrier.unique())
        carriers = append_supertype_carriers(carriers)
        carriers = remove_excluded_carriers(carriers)
        lookup[c] = carriers

    if len(n.storage_units) > 0:
        lookup["StorageUnit"] += [""]

    # lines and links
    lookup["Transmission"] = [""]
    if snakemake.config["lines_and_links_separate"]:
        lookup["Line"] = [""]
        lookup["Link"] = [""]

    return lookup


def create_mga_lookup_groups(n):

    lookup = {}
    for c in ["Generator", "StorageUnit"]:
        countries = list(np.unique([b[:2] for b in n.df(c).bus]))
        carriers = list(n.df(c).loc[n.df(c).p_nom_extendable].carrier.unique())
        carriers = append_supertype_carriers(carriers)

        combinations = []
        for country in countries:
            country_carriers = list(
                n.df(c)
                .loc[n.df(c).p_nom_extendable & [country in x for x in n.df(c).index]]
                .carrier.unique()
            )
            country_carriers = append_supertype_carriers(country_carriers)
            for carrier in country_carriers:
                combinations.append(" ".join([country, carrier]))

        lookup[c] = combinations

    lookup["Line"] = list(n.lines.apply(country_pair, axis=1).unique())
    lookup["Link"] = list(n.links.apply(country_pair, axis=1).unique())
    lookup["Transmission"] = list(
        np.unique(np.concatenate((lookup["Line"], lookup["Link"])))
    )
    if not snakemake.config["lines_and_links_separate"]:
        lookup["Line"] = []
        lookup["Link"] = []

    return lookup


def entries_from_lookup(component, elements):
    return [
        "+".join([component, e, sense]) for e in elements for sense in ["min", "max"]
    ]


if __name__ == "__main__":

    n = pypsa.Network(snakemake.input[0])

    if snakemake.wildcards.category == "hypercube":
        mga_lookup = create_mga_lookup_hypercube(n)
    elif snakemake.wildcards.category == "groups":
        mga_lookup = create_mga_lookup_groups(n)
    elif snakemake.wildcards.category == "totals":
        mga_lookup = create_mga_lookup_totals(n)
    else:
        mga_lookup = {}

    entries = []
    for component, elements in mga_lookup.items():
        entries += entries_from_lookup(component, elements)

    with open(snakemake.output[0], "w") as f:
        for entry in entries:
            f.write(entry + "\n")
