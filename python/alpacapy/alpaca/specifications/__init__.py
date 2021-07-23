""" Alpaca specifications.

The alpaca specifications provide two classes that allow to manipulate an Alpaca executable or the inputfile. New compile settings or inputfile settings
should be added in this module.
"""

# Classes and functions
from .inputfile_specifications import InputfileSpecifications
from .user_specifications import UserSpecifications
from .output_variables import OutputVariables

# Data for wildcard import (from . import *)
__all__ = [
    "InputfileSpecifications",
    "UserSpecifications",
    "OutputVariables"
]
