""" Alpaca testsuite definitions.

Different definition classes for the testsuite.
"""

# Classes and functions
from .environment import Environment
from .executable_setup import ExecutableSetup
from .folder_setup import FolderSetup
from .result_status import ResultStatus
from .testcase_information import TestcaseInfo
from .testsuite_information import TestsuiteInfo

# Data for wildcard import (from . import *)
__all__ = [
    "Environment"
    "ExecutableSetup"
    "FolderSetup"
    "ResultStatus"
    "TestcaseInfo"
    "TestsuiteInfo"
]
