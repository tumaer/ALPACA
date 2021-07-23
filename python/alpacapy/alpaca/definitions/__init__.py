""" Alpaca defintions.

Helper classes to provide a general interface for the different Alpaca specifications.
"""

# Classes and functions
from .inputfile_tag import InputfileTag
from .specifications_base import SpecificationsBase
from .tag_base import TagBase
from .user_specification_file import UserSpecificationFile
from .user_specification_tag import UserSpecificationTag
from .output_variable_tag import OutputVariableTag

# Data for wildcard import (from . import *)
__all__ = [
    "InputfileTag",
    "SpecificationsBase",
    "TagBase",
    "UserSpecificationFile",
    "UserSpecificationTag",
    "OutputVariableTag"
]
