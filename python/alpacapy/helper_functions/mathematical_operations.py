#!/usr/bin/env python3
# Python modules
from typing import List, Tuple, Dict, Union, Optional, Any, Type, IO
import numpy as np


def get_relative_error(value: float, ref_value: float) -> float:
    """ Computes the relative error between tow values.

    Parameters
    ----------
    value : float
        The value for which the relative error is computed.
    ref_value : float
        The reference value against which the value is compared.

    Returns
    -------
    float
        The relative error between both.
    """
    return (value - ref_value) / ref_value


def compute_relative_error_norms(solution: np.array, exact: np.array, cell_volume: np.array) -> np.array:
    """ Computes the relative errors norms for the L1, L2 and Linf norm between a numerical solution and its exact counter part.

    Parameters
    ----------
    solution : np.array
        The numerical solution.
    exact : np.array
        The exact solution.
    cell_volume : np.array
        The size of each numerical cell.

    Returns
    -------
    np.array
        Array with L1, L2 and Linf norm.

    Raises
    ------
    ValueError
        In case the input array are of different sizes.
    """
    if cell_volume.size != solution.size != exact.size:
        raise ValueError("Array dimensions do not match. Norms cannot be computed")
    l1 = np.sum(np.abs(solution / exact - 1.0) * cell_volume)
    l2 = np.sum((np.abs(solution / exact - 1.0) * cell_volume) ** 2.0) ** 0.5
    li = np.max(np.abs(solution / exact - 1.0) * cell_volume)
    return np.array([l1, l2, li], dtype=np.float64)


def compare_to_reference_data(computed_value: float, reference_value: float) -> float:
    """ Compares a computed value to its reference value and gives the percentage change to the reference value.

    Parameters
    ----------
    computed_value : float
        The value that is checked against reference data.
    reference_value : float
        The reference value.

    Returns
    -------
    float
        The perecentage of the change between computed and reference value (always positive value).
    """
    percentage_change = (computed_value - reference_value) / reference_value
    return percentage_change


def get_percentage(value: float) -> float:
    """ Converts a float value into its percentag value.

    Parameters
    ----------
    value : float
        The value that should be converted.

    Returns
    -------
    float
        The float in percentage, e.g. 0.5 -> 50
    """
    return round(value * 100, 2)
