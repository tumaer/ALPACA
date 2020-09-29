""" Alpaca functions.

Defines all functions that interact directly with the alpaca framework. This includes the generation of executables, run simulations or manipulate
Alpaca related files, such as log files or inputfiles.
"""

# Additional folder
from . import (
    definitions,
    specifications
)

# Classes and functions
from .create_executable import create_executable
from .obtain_runtime_information import obtain_runtime_information
from .remove_volatile_log_information import remove_volatile_log_information
from .run_alpaca import run_alpaca

# Data for wildcard import (from . import *)
__all__ = [
    "definitions",
    "specifications",
    "create_executable",
    "obtain_runtime_information",
    "remove_volatile_log_information",
    "run_alpaca"
]
