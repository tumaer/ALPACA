""" Alpaca testsuite.

Gives the possibility to run the Alpaca testsuite to check new implementations on consistency. Furthermore, it can be used to run parameter variations.
"""
# Additional folders
from . import (
    data_definitions,
    definitions,
    test_cases
)
# Classes and functions
from .testcase import Testcase
from .testsuite import Testsuite

# Data for wildcard import (from . import *)
__all__ = [
    "data_definitions"
    "definitions"
    "test_cases"
    "Testcase"
    "Testsuite"
]
