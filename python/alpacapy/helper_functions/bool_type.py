#!/usr/bin/env python3
# Python modules
from typing import List, Tuple, Dict, Union, Optional, Any, Type, IO
import numpy as np


class BoolType:
    """  Wrapper class for boolean types.

    Defines a wrapper class for boolean variables. It only allows variables to be boolean that are "true", "1", "false", "0" or any
    upper case version of those variables. It is required to specifically reject all other variables. The standard boolean class e.g.
    turns "Test" into True, which is not appropriate.

    Attributes
    ----------
    __data_type : np.dtype
        The data type of the class as numpy boolean.
    """
    __data_type = np.dtype(bool)

    @classmethod
    def __repr__(cls) -> str:
        """ Implementation of the built-in repr function """
        return str(cls.__data_type)

    @classmethod
    def __str__(cls):
        """ Implementation of the built-in str function """
        return str(cls.__data_type)

    @classmethod
    def type(cls, value: Any) -> bool:
        """ Gives the converted boolean value of a given input value.

        Parameters
        ----------
        value : Any
            The value that should be converted into bool.
        Returns
        -------
        bool
            True or False depending on the given input.
        Raises
        ------
        TypeError
            If the input type is not convertable to boolean.
        """
        # Convert the entry to string (if not possible, it already raises error)
        value = str(value).strip().lower()
        # Check if it fulfills the two criteria
        if value in ["true", "1"]:
            return True
        elif value in ["false", "0"]:
            return False
        else:
            raise TypeError("Cannot convert value into bool")
