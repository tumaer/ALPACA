#!/usr/bin/env python3
# Python modules
from typing import List, Tuple, Dict, Union, Optional, Any, Type, IO
from enum import Enum


class DataFileSuffix(Enum):
    """ The DataFileSuffix provides all possible file suffixes used for the testsuite. This allows proper file creation for all different cases.

    Parameters
    ----------
    Enum :
        Base class of the DataFileSuffix
    """
    # The file (suffix) where the relative error data is read from or written to.
    errors = "errors"
    # The file (suffix) where executable data is placed.
    executables = "executables"
    # The file (suffix) where runtime information is read from or written to.
    runtimes = "runtimes"
    # The suffix for the testcase
    testcase = "testcase"

    def __str__(self) -> str:
        """ Implementation of the built-in str function """
        return self.value

    @property
    def status(self) -> str:
        """ Adds an additional status suffix to the file suffix.

        Returns
        -------
        str
            The file suffix with _status
        """
        return self.value + "_status"

    @property
    def reference(self) -> str:
        """ Adds an additional reference suffix to the file suffix.

        Returns
        -------
        str
            The file suffix with _reference
        """
        return self.value + "_reference"
