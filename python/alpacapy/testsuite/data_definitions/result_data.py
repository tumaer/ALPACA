#!/usr/bin/env python3
# Python modules
from typing import List, Tuple, Dict, Union, Optional, Any, Type, IO
import pandas as pd
import numpy as np
import itertools
# alpacapy modules
from alpacapy.alpaca.specifications.inputfile_specifications import InputfileSpecifications
from alpacapy.testsuite.data_definitions import data_functions as data_func
from alpacapy.helper_functions import string_operations as so
from alpacapy.helper_functions import check_operations as co


class ResultData:
    """ Class to store all result and variation data of a test case.

    The ResultData class stores all data as the result from a test case. It handles writing, adding and removing of data for different
    variations of the test case. This includes varying inputfiles, executables and other test case specific variations. It further provides
    access to the data based on unique case IDs. It stores the data into a pandas.DataFrame, where the IDs are the indices and the columns describe
    different variables (e.g., inputfile). The class then defines certain functionalities to query specific data from the dataframe.

    Attributes
    ----------
    __executable_name : str
       The name of the column where the used executable name is placed.
    __inputfile_name : str
       The name of the column where the used inputfile name is placed.
    __next_variation_id : int
       The unique ID that is appended next to the data frame when new cases are added.
    __variation_keys : List[str]
       All keys of the varied variables for a test case
    __variation_type : List[Type]
       The variable type for each variation key.
    __variation_format : List[str]
       The format for each variation key, e.g. "{:16d}".
    __variation_data : pd.DataFrame
       A data frame holding a full set of variation for a single executable/inputfile combination.
    __result_data : pd.DataFrame
       The data frame holding all result data for a test case (several sets of variation data).
    """
    __executable_name = "Executable"
    __inputfile_name = "Inputfile"

    def __init__(self, variations: Dict[str, List[Any]] = {}, variations_type: Dict[str, Type] = {}, variations_format: Dict[str, str] = {},
                 index_based: bool = False) -> None:
        """ Constructor

        Parameters
        ----------
        variations : Dict[str,List[Any]], optional
            The dictionary holding for each key a list of varied values, by default {}.
        variations_type : Dict[str,Type], optional
            The dictionary holding for each variation key the type, by default {}.
        variations_format : Dict[str,str], optional
            The dictionary holding for each variation key the formatting style, by default {}.
        index_based : bool
            Flag whether the variations are done index based or not. Index based variations take the 1st element of all keys.
            Non index-based variations create all possible combinations of all list entries of the different keys.
            Example: "MaximumLevel" : [8,16], "NumberOfRanks" : [2,4]
            => Non-index base: (8,2), (8,4), (16,2), (16,4)
            => Index based: (8,2), (16,4)
        Raises
        ------
        ValueError, TypeError
            See consitency function for reference.
        """
        # Make consistency checks if variations are given
        if variations:
            self.__check_variation_consistency(variations, variations_type, variations_format, index_based)

        # Define the member variables and modify them through constructor
        self.__next_variation_id = 0
        self.__variation_keys = list(variations.keys()) + [self.__inputfile_name, self.__executable_name]
        self.__variation_type = {key: np.dtype(types) for key, types in variations_type.items()}
        self.__variation_type.update({self.__inputfile_name: np.dtype(str), self.__executable_name: np.dtype(str)})
        self.__variation_format = dict(variations_format)
        self.__variation_format.update({self.__inputfile_name: "{:s}", self.__executable_name: "{:s}"})
        # The variation data holds the block specifying all different cases for a single inputfile/executable combination
        self.__variation_data = data_func.create_dataframe_from_variations(variations, {key: "" for key in [self.__inputfile_name, self.__executable_name]},
                                                                           index_based)
        # The result data holds all cases for all different inputfile executable combinations
        self.__result_data = pd.DataFrame([], columns=self.__variation_keys)

    def __check_variation_consistency(self,
                                      variations: Dict[str, List[Any]] = {},
                                      variations_type: Dict[str, Type] = {},
                                      variations_format: Dict[str, str] = {},
                                      index_based: bool = False) -> None:
        """ Checks the consitency of the variational data the dataframe should be created of.

        Parameters
        ----------
        variations : Dict[str,List[Union[str,int]]]
            [description]
        index_based : bool, optional
            [description], by default False
        Raises
        ------
        ValueError
            1. If the keys of the variation dict do not coincide with keys of the type/format dictionaries.
            2. If the type format is not provided but the variation dictionary.
            3. If index based variation is chosen and not all variation lists have the same length.
            4. If the variations are not compliant to the inputfile tags (except NumberOfRanks).
        TypeError
            If the type or format do not coincide with the given values provided.
        """
        # Check that the names for the variations are only variables that can be modified via the InputfileSpecifications or NumberOfRanks
        co.check_list_element_existence(list(variations.keys()), list(InputfileSpecifications().keys()) + ["NumberOfRanks"],
                                        err_msg="Only inputfile specification keys or 'NumberOfRanks' are allowed in variation dict")
        # Check that the variations type and format exists for each variation
        if not variations_type or not variations_format:
            raise ValueError("The string format and type of each variation key must be specified if any is specified")
        # Check that all keys of the variations and type/format dict are the same.
        co.check_list_unique_elements(list(variations.keys()), list(variations_format.keys()),
                                      err_msg="The keys of the variations and the variations format must be the same")
        co.check_list_unique_elements(list(variations.keys()), list(variations_type.keys()),
                                      err_msg="The keys of the variations and the variations type must be the same")
        # Check the instance of the type/format dictionary
        co.check_list_element_instance(list(variations_format.values()), str, err_msg="The format dict must contain only 'str' instances")
        co.check_list_element_instance(list(variations_type.values()), type, err_msg="The type dict must contain only 'type' instances")
        # Check that all values of the variation dictionary can be converted into the provided type and format
        for (name, formats), types in zip(variations_format.items(), variations_type.values()):
            for value in variations[name]:
                if not co.check_format(value, formats):
                    raise TypeError("The specified format for " + name + " is not compatible, type of variable is " + str(type(value)))
                if not co.check_type(value, types):
                    raise TypeError("The specified type for " + name + " is not compatible, type of variable is " + str(type(value)))
        # Check if index based variations are possible
        if index_based:
            co.check_dict_for_index_based_variations(variations, err_msg="Index based variation not possible for variations. "
                                                     "Number of elements must be the same for each value of variation dict.")

    def __repr__(self) -> str:
        """ Implementation of the built-in repr function """
        string = type(self).__name__ + ":\n"
        string += "Column formats: " + str(self.__variation_format) + "\n"
        string += "Column types  : " + str(self.__variation_type) + "\n"
        string += "All data      : \n" + str(self.__result_data)
        return string

    def __getitem__(self, case_id: int) -> pd.Series:
        """ Allows accessing a single row of the result data through indexing.

        Parameters
        ----------
        case_id : int
            The unique case ID that should be accessed.
        Returns
        -------
        pd.Series
            The row of the result data as a pandas Series.
        Raises
        ------
        IndexError
            If the given case ID does not exist.
        """
        data_func.check_dataframe_index(self.__result_data, case_id)
        # Return the list (explicitly used iloc here to prohibit writing to the returned value, which would change the data inside, too)
        # iloc raise a copy error if trying to set a value to the returned series
        return self.__result_data.iloc[case_id]

    def __len__(self) -> int:
        """ Gives the number of cases present in the current data """
        return len(self.__result_data.index)

    @property
    def executable_name(self) -> str:
        """ Gives the executable name as property function """
        return self.__executable_name

    @property
    def inputfile_name(self) -> str:
        """ Gives the inputfile name as property function """
        return self.__inputfile_name

    def set_value(self, case_id: int, col: str, value: Any) -> None:
        """ Assigns a value to a certain case and column (variation).

        Parameters
        ----------
        case_id : int
            The unique case ID that should be accessed.
        col : str
            The column (variation variable) to be accessed.
        value : Any
            The value to be set.
        Raises
        ------
        IndexError
            If the case ID does not exist.
        KeyError
            If the column does not exist.
        """
        # Sanity check that the case id and column exists
        data_func.check_dataframe_index(self.__result_data, case_id)
        data_func.check_dataframe_column_names(self.__result_data, col)
        # Set the elements (type error is raised in case it does not work with the given value)
        self.__result_data.loc[case_id, col] = self.__variation_type[col].type(value)

    def itercases(self) -> Tuple[int, pd.Series]:
        """ Allows iteration through all cases row-wise.

        Returns
        -------
        Tuple[int,pd.Series]
            The case ID and row entries as pandas series.
        """
        return self.__result_data.iterrows()

    def case_ids(self, inputfile: str = None, executable: str = None) -> List[int]:
        """ Gives all case IDs for executable/inputfile combination.

        Gives all case ids where data can be found for a given inputfile and executable. If no inputfile or executable is specified,
        all indices are returned. If the inputfile and/or executable does not exist currently, a set with all variations for this combination
        is created.

        Parameters
        ----------
        inputfile : str, optional
            The inputfile for which the case ids should be returned, by default None.
        executable : str, optional
            The executable for which the case ids should be returned, by default None.
        Returns
        -------
        List[int]
            All case IDs specifying the variations for the desired data.
        Raises
        ------
        ValueError
            If no cases are found for the inputfile/executable combination and one of them is not specified. Then no set of variation can be generated.
        """
        # If no inputfile or executable is specified, return the full list of indices
        if inputfile is None and executable is None:
            return list(range(0, len(self.__result_data.index)))
        # Specify the columns that should be checked
        column_names = [] if inputfile is None else [self.__inputfile_name]
        column_names += [] if executable is None else [self.__executable_name]
        column_comparison = [] if inputfile is None else [inputfile]
        column_comparison += [] if executable is None else [executable]

        # Check if the combination of inputfile and executable already exists. If so, return the variation indices.
        indices_ref = self.__result_data.loc[(self.__result_data[column_names] == column_comparison).all(axis=1)].index
        # indices_ref = np.where( np.all( self.__result_data[ column_names ].to_numpy() == column_comparison, axis = 1 ) )
        if len(indices_ref) != 0:
            case_ids = self.__result_data.index[indices_ref]
            return case_ids.tolist()
        # If combination does not exist but only one identifier is given raise an error
        elif inputfile is None or executable is None:
            raise ValueError("A new set of variations can only be created if an inputfile AND executable is given")

        # Get the number of elements before new elements are added and the number of elements that will be added
        case_ids = list(range(self.__next_variation_id, self.__next_variation_id + len(self.__variation_data.index)))
        # Overwrite the data in the variation dataframe for the inputfile and executable in all rows
        self.__variation_data[[self.__inputfile_name, self.__executable_name]] = [inputfile, executable]
        # Set the indices to unique ids
        self.__variation_data.set_index(pd.Series(case_ids), inplace=True)
        # Add the new variations to the overall data frame (NOTE: ignore_index must be set to False to ensure that unique IDs are provided always)
        self.__result_data = pd.concat([self.__result_data, self.__variation_data], ignore_index=False)
        # Set the unique ID to the new value
        self.__next_variation_id = case_ids[-1] + 1
        # Return a list of all variation indices that have been added
        return case_ids

    def case_variation(self, case_id: int, only_variation: bool = False) -> pd.Series:
        """ Gives the variation information for a given case ID.

        Gives the variation for a given ID of the data. The variation consists of all data required to identify the case uniquely.
        This includes the inputfile, executable and test case specific variation information. The case variation information is specified during the class
        __init__(). Additional information, that is later added through add_column(), is not included in the case variation.

        Parameters
        ----------
        case_id : int
            The unique case ID that should be accessed.
        only_variation : bool, optional
            Flag whether only the variation should be included (without inputfile and executable), by default False.
        Returns
        -------
        pd.Series
            The variation data of the test case as pandas Series.

        Raises
        ------
        IndexError
            If the case ID is not found in the DataFrame.
        """
        data_func.check_dataframe_index(self.__result_data, case_id)
        keys = [key for key in self.__variation_keys if key not in ["Inputfile", "Executable"]] if only_variation else self.__variation_keys
        return self.__result_data.loc[case_id, keys]

    def case_variation_string(self, case_id: int) -> str:
        """ Gives a proper string for the variation of a single case.

        The variation string consists of the inputfile, executable and additional variation data. The output is dependent if variation data is present or not.
        1. Inputfile + Executable => "inputfile" with "executable"
        2. Inputfile + Executable + Variation => "inputfile" for variation1 = value1, variation2 = value2, ... with "executable"

        Parameters
        ----------
        case_id : int
            The unique case ID that should be accessed.
        Returns
        -------
        str
            The variation string.
        Raises
        ------
        IndexError
            If the case ID is not found in the DataFrame.
        """
        data_func.check_dataframe_index(self.__result_data, case_id)
        # Get the series for the case id
        series = self.__result_data.loc[case_id, self.__variation_keys]
        keys = [key for key in self.__variation_keys if key not in [self.__inputfile_name, self.__executable_name]]
        return series[self.__inputfile_name] + " with " + series[self.__executable_name] if series[keys].empty else \
            series[self.__inputfile_name] + " for " + ", ".join(index + " = " + str(value)
                                                                for index, value in series[keys].items()) + " with " + series[self.__executable_name]

    def add_column(self, col_name: str, col_type: type, col_format: str) -> None:
        """ Adds a colum to the result data (e.g. for additional result data).
        Parameters
        ----------
        col_name : str
            The name of the column that should be added.
        col_type : type
            The type of the column (e.g., float or int).
        col_format : str
            The format of the column that is used for writing (e.g., {:15.8f} for float).
        Raises
        ------
        ValueError
            If the type and the format do not coincide.
        Notes
        -----
        Whereas float values are not allowed for the variation data of the cases, floating values are allowed for result data
        since no comparison to those valuesis required.
        """
        # Convert the type to numpy type
        col_type = np.dtype(col_type)
        # Check the format
        if not co.check_format(col_type.type(0), col_format):
            raise ValueError("Column type and format do not coincide")
        # Add the column
        self.__result_data[col_name] = col_type.type(0)
        self.__variation_data[col_name] = col_type.type(0)
        self.__variation_type[col_name] = col_type
        self.__variation_format[col_name] = col_format

    def find_case(self, variation: pd.Series) -> Optional[List[int]]:
        """ Finds the corresponding case for a given case variation.

        Parameters
        ----------
        variation : pd.Series
            A variation series that holds the information for which the case should be found.
        Returns
        -------
        Optional[List[int]]
            The unique case ID for this user specifications (or multiple if found).
            None if variation contains data not contained in the data frame or no case is found.
        Notes
        -----
        If the variation contains values that are later added through the add_column function, it cannot be ensured that comparison works,
        if floating values are used for those columns.
        """
        # Find the case in the dataframe
        return data_func.find_case_in_dataframe(self.__result_data, variation)

    def del_column(self, *col_names: Tuple[str]) -> None:
        """ Deletes a single or multiple columns from the dataframe.

        Parameters
        ----------
        col_names : Tuple[str]
           The unpacked list of column names that should be removed (lists must be unpacked in function call, otherwise error).
        Notes
        -----
        If columns, that specify variation data, are removed, the class could throw errors at some other places (e.g., find case.)
        """
        # Only deletes column if present, otherwise do nothing is done
        for col_name in col_names:
            if col_name in self.__result_data:
                del self.__result_data[col_name]
                del self.__variation_data[col_name]

    def del_case(self, *case_ids: Tuple[int]) -> None:
        """ Deletes a single or multiple cases from the dataframe.

        Parameters
        ----------
        case_ids : Tuple[int]
           The unpacked list of case IDs that should be removed (lists must be unpacked in function call, otherwise error).
        """
        # Only deletes case if present, otherwise nothing is done
        for case_id in case_ids:
            if case_id in self.__result_data.index:
                self.__result_data.drop(case_id, inplace=True)

    def to_csv(self, name: str, *columns: Tuple[str]) -> IO:
        """ Writes the current executable data to a csv file.

        During the creation of the csv file, all columns (including header) are formatted based on the maximum length of all values per column.

        Parameters
        ----------
        name : str
            The ABSOLUTE path to the csv file where the data is written to.
        columns : str
            The unpacked list of column name(s) that are written to the file. The variation identifiers are always written.
            If no columns are given all data is written.
        Raises
        ------
        IOError
           If columns are specified that are not in the dataframe (from pandas).
        """
        # Check if columns exists
        for column in columns:
            if column not in self.__result_data:
                raise IOError("Cannot write column '" + column + "' to csv-file. Column does not exist.")
        # Only take columns that are not part of the variations definitions
        columns = [column for column in columns if column not in self.__variation_keys]
        # Determine all columns that should be written
        columns_to_be_written = columns + self.__variation_keys if columns else list(self.__result_data.columns)
        # First copy the DataFrame to a temporary object to allow formatting
        dataframe_to_write = self.__result_data[columns_to_be_written].copy()
        # Compute the maximum length of inputfiles and executables and add this format to the variation_format. The +5 is a simple offset to add whitespace
        # before the entry to improve file readability.
        for column in ["Inputfile", "Executable"]:
            self.__variation_format[column] = data_func.get_string_column_max_length_format(dataframe_to_write, column)
        # Format the dataframe
        dataframe_to_write = data_func.format_dataframe(dataframe_to_write, self.__variation_format, include_header=True)
        # Write the dataframe to file
        dataframe_to_write.to_csv(name, sep=",", header=True, index=False, encoding='utf-8', na_rep="NaN")
