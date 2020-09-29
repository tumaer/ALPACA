#!/usr/bin/env python3
# Python modules
from typing import List, Tuple, Dict, Union, Optional, Any, Type, IO
import numpy as np
# alpacapy modules
from alpacapy.helper_functions.bool_type import BoolType


def bool_to_active(value: bool) -> str:
    """ Gives Active or Deactivated depending on the boolean value.

    Parameters
    ----------
    value : bool
        The boolean value that should be converted.
    Returns
    -------
    str
        Active if True, otherwise Deactivated.
    """
    return "Active" if value else "Deactivated"


def string_to_bool(value: str) -> bool:
    """ Converts a string into a boolean value only for specific conditions.
    In contrast to built-in string-to-bool conversion where every set string becomes True.

    Parameters
    ----------
    value : str
        The value that should be converted.
    Returns
    -------
    bool
        True or False.
    Raises
    ------
    TypeError
        If the value is not convertable to boolean.
    """
    return BoolType.type(value)


def string_to_float(value: str, default: float) -> float:
    """ Converts a string into float value.

    Parameters
    ----------
    value : str
        The value that should be converted.
    default : float
        Fallback default value if string not convertible.

    Returns
    -------
    float
        The converted/default value.
    """
    try:
        return float(value)
    except ValueError:
        return default


def dim_to_str(dim: int, xml: bool = False) -> str:
    """ Converts a dimension into a proper string.

    Parameters
    ----------
    dim : int
        The dimension to be converted.
    xml : bool, optional
        Flag whether xml style is used, by default False.
    Returns
    -------
    str
        The converted string.
    """
    if dim == 1:
        return "one" if xml else "One"
    elif dim == 2:
        return "two" if xml else "Two"
    else:  # dim ==3
        return "three" if xml else "Three"


def convert_to_percentage(value: Optional[float] = None, width: Optional[int] = None, precision: int = 0) -> str:
    """ Converts a floating value to a percentage string.

    Parameters
    ----------
    value : Optional[float], optional
        The value to be converted, by default None.
    width : Optional[int], optional
        The width of the string generated, by default None.
    precision : int, optional
        The precision of the percentag value, by default 0.
    Returns
    -------
    str
        The percentags string.
    """
    if value is not None and not isinstance(value, str):
        if width is None:
            width = ""
        return ("{:" + str(width) + "." + str(precision) + "f}").format(round(value * 100, 2))
    else:
        return value


def convert_value_to_string(value: Any, width: Optional[int] = None, precision: int = 8) -> str:
    """ Converts a general value into a proper string.

    Parameters
    ----------
    value : Any
        The value to be converted.
    width : Optional[int], optional
        The width of the string, by default None.
    precision : int, optional
        The precision of floating values, by default 8.
    Returns
    -------
    str
        The converted value.
    """
    # Differ between float, int and str
    if isinstance(value, float):
        width = "" if width is None else width
        return ("{:" + str(width) + "." + str(precision) + "e}").format(value)
    elif isinstance(value, int):
        return str(value) if width is None else str(value).rjust(width)
    else:
        return value if width is None else value.rjust(width)


def cut_string(string: str, width: int) -> List[str]:
    """ Cuts the string into substrings of a maximum width.

    Parameters
    ----------
    string : str
        Thr string that should be cut.
    width : int
        The width at which the string is cut.
    Returns
    -------
    List[str]
        The list with all substrings
    """
    return [string[i:i + width] for i in range(0, len(string), width)]


def remove_vowels(string: str) -> str:
    """ Removes all vowels from a string.

    Parameters
    ----------
    string : str
        The string where the vowels should be removed.
    Returns
    -------
    str
        The modified string.
    """
    vowels = ["a", "b", "c", "d", "e"]
    return "".join([char for char in string if char.lower() not in vowels])
