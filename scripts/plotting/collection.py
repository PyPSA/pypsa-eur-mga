"""
Collects various plotting functions.
"""

__author__ = "Fabian Neumann (KIT)"
__copyright__ = "Copyright 2019-2020 Fabian Neumann (KIT), GNU GPL 3"

from .bar import plot_capacity_bar, plot_cost_bar
from .boundaries import plot_space, plot_space_presentation
from .boxplots import plot_curated_boxplots, plot_plain_boxplot
from .dominance import plot_dominance
from .gini import plot_lorentz, plot_gini
from .map import plot_network
from .multbar import plot_bar_collection
from .pie import plot_energy_pie
from .violins import plot_violins
from .correlations import (
    plot_capacity_correlation,
    plot_energy_correlation,
    plot_correlation,
)
