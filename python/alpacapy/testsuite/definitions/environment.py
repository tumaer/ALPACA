#!/usr/bin/env python3
# Python modules
from typing import List, Tuple, Dict, Union, Optional, Any, Type, IO
from enum import Enum


class Environment(Enum):
    """ The Environment provides all possible choices of the environment, where the testsuite can be run.

    Parameters
    ----------
    Enum :
        Base class of the Environment.
    """
    merge = {"Name": "MERGE", "CompileCores": 16, "RunRanks": [5, 13, 28]}
    aer = {"Name": "AER", "CompileCores": 8, "RunRanks": [4, 7, 12]}
    sumuc = {"Name": "SUMUC", "CompileCores": 32, "RunRanks": [17, 35, 48]}
    cluster = {"Name": "CLUSTER", "CompileCores": 10, "RunRanks": [16, 28, 28]}
    custom = {"Name": "CUSTOM", "CompileCores": 1, "RunRanks": [2, 2, 2]}
    not_known = {"Name": "UNKNOWN", "CompileCores": 0, "RunRanks": [1, 1, 1]}

    def __str__(self) -> str:
        """ Implementation of the built-in str function. Gives the environment name. """
        return str(self.value["Name"][0].upper() + self.value["Name"].lower()[1:])

    @classmethod
    def get_env(cls, name: str) -> 'Environment':
        """ Returns the correct environment for a given name.

        Parameters
        ----------
        name : str
            The name of the environment for which the correct enum variable should be obtained.
        Returns
        -------
        Environment
            The correct environment enumeration variable.
        """
        # Convert the name to upper case letters
        name = name.upper()
        if name in list(map(lambda c: c.value["Name"], cls)):
            return cls[name.lower()]
        else:
            return cls.not_known

    @classmethod
    def get_custom_env(cls, number_of_ranks: List[int]) -> 'Environment':
        """ Creates the custom environment with the desired number of ranks.

        Parameters
        ----------
        number_of_ranks : List[int]
            The number of ranks used for the custom environment.
        """
        # Check if all number of ranks are given
        if len(number_of_ranks) != 3:
            ValueError("The number of ranks for the custom environment must be 3")
        cls.custom.value["RunRanks"] = number_of_ranks
        cls.custom.value["CompileCores"] = max(number_of_ranks)
        return cls.custom

    def compile_cores(self) -> int:
        """ Gives the number of compile cores.

        Returns
        -------
        int
            The number of compile cores.
        """
        return self.value["CompileCores"]

    def number_of_ranks(self, dim: Optional[int] = None) -> Union[List[int], int]:
        """ Gives the number of ranks to run a simulation.

        Parameters
        ----------
        dim : Optional[int], optional
            dimension for which the number of ranks should be given, by default None. If None all are given.

        Returns
        -------
        Union[List[int], int]
            All ranks or the dimensional rank.
        """
        return self.value["RunRanks"][dim - 1] if dim is not None else self.value["RunRanks"]
