# Python modules
from typing import List, Tuple, Dict, Union, Optional, Any, Type, IO
import pandas as pd
import itertools


def check_dataframe_index(dataframe: pd.DataFrame, index: int) -> None:
    """ Checks if a certain index is present in the dataframe.

    Parameters
    ----------
    dataframe : pd.DataFrame
        The dataframe to be checked.
    index : int
        The index to be checked.
    Raises
    ------
    IndexError
        If the index is not present in the dataframe.
    """
    if index not in dataframe.index:
        raise IndexError("The given index " + str(case_id) + " does not exist in dataframe")


def check_dataframe_column_names(dataframe: pd.DataFrame, *columns: str) -> None:
    """ Checks if the column(s) are present in the dataframe.

    Parameters
    ----------
    dataframe : pd.DataFrame
        The dataframe to be checked.
    columns : str
        The column(s) to be checked (unpacked list).
    Raises
    ------
    KeyError
        If the column is not present.
    """
    if any(column not in dataframe for column in columns):
        raise ValueError("The column '" + column + "' cannot be found in the dataframe")


def create_dataframe_from_variations(variation_dict: Dict[str, List[Any]], default_dict: Dict[str, List[Any]] = {}, index_based: bool = False) -> pd.DataFrame:
    """ Creates a pandas dataframe based on a given dictionary of variations.

    Parameters
    ----------
    variation_dict : Dict[str,List[Any]]
        The dictionary holding all values that should be varied.
    default_dict : Dict[str,List[Any]], optional
        The dictionary adding default values to all combinations, by default {}.
    index_based : bool, optional
        Flag whether index based variations should be created, by default False.
    Returns
    -------
    pd.DataFrame
        The fully created pandas dataframe
    """
    if index_based:
        # Loop through all keys and take the n-th element and combine them into a tuple
        combinations = [tuple(value[i] for value in variation_dict.values()) + tuple(default_dict.values())
                        for i in range(0, len(list(variation_dict.values())[0]))]
    else:
        # Take all lists of the keys, add an empty string for the inputfile and executable and generate the combinations
        combinations = [variation_list for variation_list in variation_dict.values()]
        # NOTE: The order of each tuple created is equal to the order of the keys
        combinations = list(itertools.product(*combinations))
        # Add the default values to the combination
        combinations = [combination + tuple(default_dict.values()) for combination in combinations]

    return pd.DataFrame(combinations, columns=list(variation_dict.keys()) + list(default_dict.keys()))


def find_case_in_dataframe(dataframe: pd.DataFrame, series: pd.Series) -> Optional[List[int]]:
    """ Finds a series in a dataframe.

    Parameters
    ----------
    dataframe : pd.DataFrame
        The dataframe where the series should be found.
    series : pd.Series
        The series that should be found.
    Returns
    -------
    Optional[List[int]]
        List with all indices that match the series. None if no is found.
    """
    # If any index in the series cannot be found in the executable data return None
    if any(index not in dataframe for index in series.index):
        return None
    # Find the indices (multiple indices are possible)
    indices = dataframe.loc[(dataframe[series.index] == series).all(axis=1)].index
    # Depending on the output return index or None
    return dataframe.index[indices].tolist() if len(indices) != 0 else None


def get_string_column_max_length_format(dataframe: pd.DataFrame, column: str) -> str:
    """ Gets the format of dataframe column of type string that gives the maximum entry length.

    Parameters
    ----------
    dataframe : pd.DataFrame
        The dataframe where column is found.
    columns : str
        The column where the format should be given
    Returns
    -------
        The string that can be used for formatting (e.g., {:>100s}). Always right adjusted with 5 white-space offset.
    """
    return "{:>" + str(max(max(dataframe[column].apply(len)), len(column)) + 5) + "s}"


def format_dataframe(dataframe: pd.DataFrame, format_dict: Dict[str, str] = None, format_type: Type = None, include_header: bool = True) -> pd.DataFrame:
    """ Formats the header of a dataframe based on the maximum length of the column entries.

    Parameters
    ----------
    dataframe : pd.DataFrame
        The dataframe that is modified. In-place modification.
    format_dict : Dict[str,str], optional
        The format dictionary that is applied on the dataframe (all columns of the dataframe must be contained), by default None.
    format_type : Type, optional
        The format type that is used, by default None.
    include_header : bool, optional
        Flag whether the header should be included.
    Returns
    -------
    pd.DataFrame
        The modifid dataframe.
    """
    if format_dict is None and format_type is None:
        raise ValueError("Either a format dict or type must be specified")
    # Format the columns based on dict or type (dict has priority)
    dataframe = dataframe.apply(lambda x: x.apply(format_dict[x.name].format)) if format_dict is not None \
        else dataframe.apply(lambda x: x.apply(format_type))
    # format the header based on the first column entry len
    if include_header:
        formatted_header = {name: ("{:>" + str(len(dataframe[name].iloc[0])) + "s}").format(name) for name in dataframe}
        dataframe.rename(columns=formatted_header, inplace=True)
    return dataframe
