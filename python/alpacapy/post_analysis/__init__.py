""" Alpaca post-analysis module

Provides modules that compare the result of Alpaca simulations against analytical or other data.
"""

# Classes and function
from .check_symmetry import check_symmetry
from .couette_flow_two_interfaces import CouetteFlowTwoInterfaces
from .oscillating_drop import OscillatingDrop
from .shear_drop_deformation import ShearDropDeformation
from .sod_analysis import SodAnalysis

# Data for wildcard import (from . import *)
__all__ = [
    "check_symmetry",
    "CouetteFlowTwoInterfaces",
    "OscillatingDrop",
    "ShearDropDeformation",
    "SodAnalysis"
]
