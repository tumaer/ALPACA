#!/usr/bin/env python3
# Python modules
from typing import List, Tuple, Dict, Union, Optional, Any, Type, IO


def check_format(value: Any, value_format: str) -> bool:
    """ Checks that a given formatting string is possible for the value.

    Parameters
    ----------
    value : Any
        The value for which the format should be checked.
    value_format : str
        The format applied on the value, e.g. "{:13.7f}".
    Returns
    -------
    bool
        True if formatting is possible, False otherwise.
    """
    try:
        value_format.format(value)
        return True
    except ValueError:
        return False


def check_type(value: Any, type_style: Type) -> bool:
    """ Checks whether a value can be converted into a given typ.

    Parameters
    ----------
    value : Any
        The value for which the type should be checked.
    type_style : Type
        The type function applied on the value.
    Returns
    -------
    bool
        True if possible, False otherwise.
    """
    try:
        type_style(value)
        return True
    except TypeError:
        return False


def check_list_element_existence(list_to_check: List[Any], ref_list: List[Any], err_msg: str = "") -> None:
    """ Checks if all elements from a list are given in a reference list.

    Parameters
    ----------
    list_to_check : List[Any]
        The list that is checked.
    ref_list : List[Any]
        The reference list.
    err_msg : str, optional
        The error message that is printed, by default "".
    Raises
    ------
    ValueError
        If any key of the list are not found in the reference list.
    """
    if any(key not in ref_list for key in list_to_check):
        raise ValueError(err_msg)


def check_list_unique_elements(list_to_check: List[Any], ref_list: List[Any], err_msg: str = "") -> None:
    """ Checks if all elements from one list a provided in another list.

    Parameters
    ----------
    list_to_check : List[Any]
        The list that is checked.
    ref_list : List[Any]
        The reference list.
    err_msg : str, optional
        The error message that is printed, by default "".
    Raises
    ------
    ValueError
        If both list have different unique elements.
    """
    if set(list_to_check) != set(ref_list):
        raise ValueError(err_msg)


def check_list_element_instance(list_to_check: List[Any], *value_types: Type, err_msg: str = "") -> None:
    """ Checks if all values in a dictionary have a certain instance type.

    Parameters
    ----------
    list_to_check : List[Any]
        The list that is checked.
    value_types : Type
        The type each list value should have. Unpacked list for multiple types.
    err_msg : str, optional
        The error message that is printed, by default "".
    Raises
    ------
    TypeError
        If any value has not the correct instance.
    """
    # Check each type individually. If any fulfills criterion, the check passes
    check_failed = True
    for value_type in value_types:
        if all(isinstance(value, value_type) for value in list_to_check):
            check_failed = False

    if check_failed:
        raise TypeError(err_msg)


def check_dict_for_index_based_variations(dict_to_check: Dict, err_msg: str = "") -> None:
    """ Checks whether a dictionary can be used to create index based variations.

    This function can only dictionary, where the values are iterables that provide the len() function.

    Parameters
    ----------
    dict_to_check : Dict
        The dictionary that is checked.
    err_msg : str, optional
        The error message that is printed, by default "".
    Raises
    ------
    TypeError
        If the values of the dict are not of type list or tuple.
    ValueError
        If the index based variation is not possible.
    """
    # First check that each key provides a list or tuple of elements (to allow iterations).
    check_list_element_instance(list(dict_to_check.values()), list, tuple, err_msg=err_msg)
    # Check that each value has the same number of elements.
    first_key = list(dict_to_check.keys())[0]
    ref_length = len(dict_to_check[first_key])
    if any(len(value) != ref_length for value in dict_to_check.values()):
        raise ValueError(err_message)
