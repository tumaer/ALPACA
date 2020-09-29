#!/usr/bin/env python3
# Python modules
from typing import List, Tuple, Dict, Union, Optional, Any, Type, IO
import pandas as pd
# alpacapy modules
from alpacapy.testsuite.data_definitions import data_functions as data_func
from alpacapy.helper_functions import check_operations as co


class ResultReferenceData:
    """ Class to store all data for reference values of a test case.

    The ResultReferenceData class stores all information about the reference values for a single result data. It holds the reference values and
    the naming of the reference values for proper accessing.

    Attributes
    ----------
    __reference_variables : List[str]
       All variables specifying reference value for a test case (e.g., DensityError).
    __reference_data : pd.DataFrame
       The data frame holding all reference data for a test case.
    """

    def __init__(self, reference_variables: List[str] = []) -> None:
        """ Constructor.

        Parameters
        ----------
        reference_variables : List[str], optional
            The variables given reference values for the test case, by default [].
        Raises
        ------
        TypeError
            If the reference variables are not strings.
        """
        # Make consistency check
        co.check_list_element_instance(reference_variables, str, err_msg="All reference value specifier must be of type str")
        # Data passed (store the values and create empty DataFrame for the reference values)
        self.__reference_variables = reference_variables
        self.__reference_data = pd.DataFrame([])

    def __repr__(self) -> str:
        """ Implementation of the built-in repr function """
        string = "Reference variables: " + str(self.__reference_variables)
        string += "\nReference data    :\n" + str(self.__reference_data)
        return string

    def __getitem__(self, case_id: int) -> pd.Series:
        """ Allows accessing a single row of the reference data through indexing.

        Parameters
        ----------
        case_id : int
            The unique case ID that should be accessed.
        Returns
        -------
        pd.Series
            The row of the reference data as a pandas Series.
        Raises
        ------
        IndexError
            1. If the given case ID does not exist.
            2. If the reference data is empty.
        """
        # Sanity check that the case id exists
        if self.__reference_data.empty:
            raise IndexError("No reference data exist to read data from. Read the reference values before.")
        data_func.check_dataframe_index(self.__reference_data, case_id)
        # Return the reference values for the given case id
        return self.__reference_data.loc[case_id, self.__reference_variables]

    @property
    def names(self) -> List[str]:
        """ allows .names of all reference variables that are used

        Returns
        -------
        List[str]
           List with all variables.
        """
        return self.__reference_variables

    @property
    def active(self) -> bool:
        """ allows .active call to check whether reference variales exists or not

        Returns
        -------
        bool
           True if any exists, False otherwise.
        """
        return True if self.__reference_variables else False

    def read_reference_values(self, file_path: str) -> None:
        """ Reads the reference values from a given file.

        Parameters
        ----------
        file_path : str
            The ABOSLUTE path to the file where the reference values can be found.
        Raises
        ------
        ValueError
            If any variable cannot be found in the reference data file.
        """
        # Read data and check if all columns of the reference variables are present
        if self.__reference_data.empty:
            self.__reference_data = pd.read_csv(file_path, sep=' *, *', skipinitialspace=True, engine="python")
            data_func.check_dataframe_column_names(self.__reference_data, *self.__reference_variables)

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
        Raises
        ------
        ValueError
            If the reference data is empty.
        """
        # Raise error if not reference data is present and if variation columns do not exists in data
        if self.__reference_data.empty:
            raise ValueError("No reference data exist to read data from. Read the reference values before.")
        # Find the case in the dataframe
        return data_func.find_case_in_dataframe(self.__reference_data, variation)
