""" alpacapy module

The alpacapy module provides all functionality to manipulate Alpaca result data, create executables, run simulations and test the framework.

Formatting according to PEP8 can be ensured using:
   pip install autopep8
   pycodestyle  --max-line-length=160 ./python (for listing format problems)
   autopep8 --in-place --aggressive --max-line-length=160 <file-to-be-modified) ( for auto-formatting, manual fixes possible)

A documentation of the module can be generated calling:
   pip install pdoc3
   pdoc --html --output-dir doc <path-to-alpacapy>
To exclude single files or folders from the documentation (e.g., if special modules are required that need to be installed manually) add those to the __pdoc__
dictionary below. Example:
    __pdoc__["alpaca"] = False (excludes the whole alpaca folder)
    __pdoc__["logger"] = False (excludes only the logger file)
"""
__pdoc__ = {}

__version__ = "0.0.1"

##
# Modules that can be loaded directly without invoking the subfolders
# Usage example: from alpacapy.logger import Logger

# The logger used for proper logging of all modules
from alpacapy.logger import Logger
from alpacapy.name_style import NameStyle

# Alpaca modules
from alpacapy.alpaca.specifications.user_specifications import UserSpecifications
from alpacapy.alpaca.specifications.inputfile_specifications import InputfileSpecifications
from alpacapy.alpaca.create_executable import create_executable
from alpacapy.alpaca.obtain_runtime_information import obtain_runtime_information
from alpacapy.alpaca.remove_volatile_log_information import remove_volatile_log_information
from alpacapy.alpaca.run_alpaca import run_alpaca

# PostProcessing
from alpacapy.post_processing.trim_output import trim_files

# PostAnalysis
from alpacapy.post_analysis.couette_flow_two_interfaces import CouetteFlowTwoInterfaces
from alpacapy.post_analysis.sod_analysis import SodAnalysis
from alpacapy.post_analysis.shear_drop_deformation import ShearDropDeformation
from alpacapy.post_analysis.check_symmetry import check_symmetry

# The full Testsuite
from alpacapy.testsuite.testsuite import Testsuite

__all__ = [
    "Logger"
    "NameStyle"
    "alpaca"
    "helper_functions"
    "post_processing"
    "post_analysis"
    "testsuite"
]
