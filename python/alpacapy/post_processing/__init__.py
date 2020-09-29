""" Alpaca post-processing module

Allows to process result data from an Alpaca simulation for further use, such as field manipulation and data extraction. It does not compare any data
to existing results, but creates new h5 files, such as trimmed data, or images, such as plots.
"""

# Classes and function
from .trim_output import trim_files

# Data for wildcard import (from . import *)
__all__ = [
    "trim_files"
]
