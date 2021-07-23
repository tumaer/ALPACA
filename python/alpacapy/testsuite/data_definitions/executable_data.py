#!/usr/bin/env python3
# Python modules
from typing import List, Tuple, Dict, Union, Optional, Any, Type, IO
import pandas as pd
import itertools
import xml.etree.ElementTree as et
# alpacapy modules
from alpacapy.logger import Logger
from alpacapy.alpaca.specifications.user_specifications import UserSpecifications
from alpacapy.helper_functions import check_operations as co
from alpacapy.helper_functions import xml_operations as xo
from alpacapy.helper_functions import string_operations as so
from alpacapy.name_style import NameStyle
from alpacapy.testsuite.data_definitions import data_functions as data_func


class ExecutableData:
    """ Class holding all information about the Testsuite executables.

    The ExecutableData class stores all data that is provided for an executable that is created during a testsuite run. It handles the writing, adding and
    removing of data for different variations of the executable, especially varying user specifications. It further provides access of the data based on
    unique variation/case IDs. Executables that are not created through the Testsuite are not handled.

    Attributes
    ----------
    __executable_data : pd.DataFrame
       A pandas DataFrame holding all information about the status of user specifications for a given executable where the column name equals the
       UserSpecification name
    """

    def __init__(self, executable_data: pd.DataFrame):
        """ Constructor

        Parameters
        ----------
        executable_data : pd.DataFrame
            The DataFrame that is stored as executable data.
        Notes
        -----
        Should not be called directly since no sanity check is done on data consistency. Use from_variations and from_csv instead.
        """
        # Copy the executable data to the instance variables
        self.__executable_data = executable_data.copy()

    @classmethod
    def __check_variation_consistency(cls, variations: Dict[str, List[Union[str, int]]], index_based: bool) -> None:
        """ Checks the consistency of the data in a variation dict.

        Parameters
        ----------
        variations : Dict[str,List[Union[str,int]]]
            The dictionary containing the variation information.
        index_based : bool, optional
            Flag whether index based variations are desired or not.
        Raises
        ------
        ValueError
            1. If the variation keys are not present in the user specifications.
            2. If index based variation is chosen and not all variation lists have the same length.
        TypeError
            If the type of the list elements does not coincide with the required data.
        """
        # Only variable names of the user specifications are allowed to be varied.
        co.check_list_element_existence(list(variations.keys()), list(UserSpecifications().keys()),
                                        err_msg="Only user specification keys are allowed in variation dict")
        # Check that the values are of type list.
        co.check_list_element_instance(list(variations.values()), list, err_msg="Value of each specification variation must be of type 'list'")
        # Check if index based variations are possible
        co.check_dict_for_index_based_variations(variations, err_msg="Index based variation not possible for user specification."
                                                 " Number of elements must be the same for each value of variation dict.")

    @classmethod
    def __create_dataframe_from_variations(cls, variations: Dict[str, List[Union[str, int]]],
                                           dim: int, alpaca_base_path: str, index_based: bool) -> pd.DataFrame:
        """ Creates the pandas dataframe for the set of variations.

        Parameters
        ----------
        variations : Dict[str,List[Union[str,int]]]
            The dictionary containing the variation information.
        dim : int
            The dimension of the variation data.
        alpaca_base_path : str
            The base path to the alpaca src folder.
        index_based : bool
            Flag whether index based variations are desired or not.
        Returns
        -------
        pd.DataFrame
            The fully created pandas dataframe.
        """
        # Read all user specifications from the file
        user_specifications = UserSpecifications()
        user_specifications.read_specifications(alpaca_base_path)
        # Create a dictionary with the default values that are not in the variation dictionary
        default_dict = {name: specification.value for name, specification in user_specifications.items() if name not in variations.keys()}
        # Convert all values from the variations to the correct values of the user specifications
        variation_dict = {key: [user_specifications[key].type(value)] for key, values in variations.items() for value in values}
        # Create the dataframe based on the given variations and default values
        executable_data = data_func.create_dataframe_from_variations(variations, default_dict, index_based)
        # All values of the user specifications have been added correctly. Now create the executable names for each variation and add it to the dataframe
        user_specifications.reset()
        executable_names = ["" for _ in range(0, len(executable_data.index))]
        for index, row in executable_data.iterrows():
            for key in variations.keys():
                user_specifications[key].value = row[key]
            executable_name = user_specifications.add_tags("ALPACA", dim, use_abbreviation=False)
            # There is a possibility that the name is to long to be used for file writing (the maximum allowed length is 255 characters
            # for most operating systems). In case the name is too long remove the vowels. If this is not sufficient throw an error.
            executable_names[index] = so.remove_vowels(executable_name) if len(executable_name) >= 255 else executable_name
        executable_data["Executable"] = executable_names

        return executable_data

    @classmethod
    def from_variations(cls, alpaca_base_path: str, dim: int,
                        user_spec_variations: Dict[str, List[Union[str, int]]], index_based: bool = False) -> 'ExecutableData':
        """ Overloaded constructor.

        This overloaded constructor call allows the instantiation of the class from a variation dictionary. Call ExecutableData.from_variations().

        Parameters
        ----------
        alpaca_base_path : str
            The relative or absolute path to the alpaca base folder, where the src folder lies.
        dim : int
            The dimension this executable setup is built of.
        user_spec_variations : Dict[str,List[Union[str,int]]]
            A dictionary built from the user specification name and a list of values that are varied.
        index_based : bool
            Flag whether the variations are done index based or not. Index-based variations take the n-th element of all keys and creates a variaton.
            Non index-based variations create all possible combinations of all list entries of the different keys.
            Example: "InternalCells" : [8,16], "HaloSize" : [2,4]
            => Non-index base: (8,2), (8,4), (16,2), (16,4)
            => Index based: (8,2), (16,4)
        Returns
        -------
        ExecutableData
           The instantiated ExecutableData class.
        Raises
        ------
        ValueError
            1. If the keys of the variation dict do not coincide with the user specification names
            2. If index based variation is chosen and not all list of specifications have the same length.
        TypeError
            If the values of the variation dict are not a list and each list entry is not a single element
        """
        # Check the consistency of the variation dict and obtain the dictionary with correctly converted values to be used for the user specifications
        cls.__check_variation_consistency(user_spec_variations, index_based)
        # All checks passed. Now create all combinations for the variations.
        return cls(cls.__create_dataframe_from_variations(user_spec_variations, dim, alpaca_base_path, index_based))

    @classmethod
    def from_csv(cls, path_to_file: str) -> 'ExecutableData':
        """ Overloaded constructor.

        This overloaded constructor call allows the instantiation of the class from an existing csv file. Call ExecutableData.from_csv().

        Parameters
        ----------
        path_to_file : str
            The ABSOLUTE path to the csv file that is read.
        Returns
        -------
        ExecutableData
           The instantiated ExecutableData class.
        """
        logger = Logger()
        # Instantiate the user specifications
        user_specifications = UserSpecifications()
        # Read the full file
        executable_data = pd.read_csv(path_to_file, sep=' *, *', skipinitialspace=True, engine="python")
        # Check that the columns coincide with the user specifications (all of them must be present) and that Executable is given
        try:
            co.check_list_element_existence(list(user_specifications.keys()) + ["Executable"], executable_data.columns.tolist())
        except ValueError:
            logger.write("The provided executable data file does not contain all user specifications and 'Executable'", color="y")
            return cls(pd.DataFrame([]))
        # Convert all columns to its appropriate type specified at the user specifications
        executable_data = executable_data.apply(lambda x: x.apply(user_specifications[x.name].type) if x.name != "Executable" else x)
        # Create the class instance
        return cls(executable_data)

    def __repr__(self) -> str:
        """ Implementation of the built-in repr function """
        string = type(self).__name__ + ":\n"
        string += "All data: \n" + str(self.__executable_data)
        return string

    def __getitem__(self, case_id: int) -> pd.Series:
        """ Allows accessing a single row of the executable data through indexing.

        Parameters
        ----------
        case_id : int
            The unique case ID that should be accessed.
        Returns
        -------
        pd.Series
            The row of the executable data as a pandas Series.
        Raises
        ------
        IndexError
            If the given case ID does not exist.
        """
        data_func.check_dataframe_index(self.__executable_data, case_id)
        # Return the list (explicitly used iloc here to prohibit writing to the returned value, which would change the data inside, too)
        # iloc raise a copy error if trying to set a value to the returned series
        return self.__executable_data.iloc[case_id]

    def __len__(self) -> int:
        """ Gives the number of cases present in the current data """
        return len(self.__executable_data.index)

    @property
    def empty(self) -> bool:
        """ Checks if the data frame is empty or not.

        Returns
        -------
        bool
            True if empty, False otherwise.
        """
        return self.__executable_data.empty

    def case_ids(self) -> List[int]:
        """ Gives all case IDs present in the current data frame.

        Returns
        -------
        List[int]
            A list with all unique case IDs.
        """
        return list(self.__executable_data.index)

    def set_value(self, case_id: int, col: str, value: Any) -> None:
        """ Assigns a value to a certain case and colum (specification).

        Parameters
        ----------
        case_id : int
            The unique case ID that should be accessed.
        col : str
            The column (user specification) to be accessed.
        value : Any
            The value to be set.
        Raises
        ------
        IndexError
            If the case ID does not exist.
        KeyError
            If the column does not exist.
        """
        # Sanity check that the case id exists
        data_func.check_dataframe_index(self.__executable_data, case_id)
        data_func.check_dataframe_column_names(self.__executable_data, col)
        # Set the elements (type error is raised in case it does not work with the given value)
        self.__executable_data.loc[case_id, col] = value

    def get_executable_name(self, case_id: int) -> str:
        """ Gives the executable name for a given case ID.

        Parameters
        ----------
        case_id : int
            The unique case ID that should be accessed.
        Returns
        -------
        str
            The executable name.
        Raises
        ------
        IndexError
            If the case ID does not exist.
        """
        data_func.check_dataframe_index(self.__executable_data, case_id)
        # Return the user specifications
        return self.__executable_data.loc[case_id]["Executable"]

    def get_user_specifications(self, case_id: int) -> UserSpecifications:
        """ Gives the user specifications for a given case ID.

        Parameters
        ----------
        case_id : int
            The unique case ID that should be accessed.
        Returns
        -------
        UserSpecifications
            The fully filled user specifications.
        Raises
        ------
        IndexError
            If the case ID does not exist.
        """
        data_func.check_dataframe_index(self.__executable_data, case_id)
        # Obtain the series for the id and return the created user specification
        return UserSpecifications.from_pandas(self.__executable_data.loc[case_id])

    def find_case(self, user_specifications: UserSpecifications) -> Optional[List[int]]:
        """ Finds the corresponding case for a given set of user specifications.

        Parameters
        ----------
        user_specifications : UserSpecifications
            The fully initialized user specifications class for which the case ID should be found. If 'None' values are present None is returned.
        Returns
        -------
        Optional[List[int]]
            The unique case ID for this user specifications (or multiple if found).
            None if the user specification class is not fully initialized or no case is found.
        """
        # Create a pandas series with the user specification data
        series = pd.Series([spec.value for spec in user_specifications.values()], index=user_specifications.keys())
        # Find the case in the dataframe
        return data_func.find_case_in_dataframe(self.__executable_data, series)

    def find_executable(self, executable_name: str) -> Optional[int]:
        """ Finds a case based on a given executable name.

        Parameters
        ----------
        executable_name : str
           The executable for which the case should be found.
        Returns
        -------
        Optional[int]
            The unique case ID. None if no is found.
        """
        index = self.__executable_data.loc[self.__executable_data["Executable"] == executable_name].index
        return None if len(index) == 0 else index[0]

    def to_csv(self, name: str) -> IO:
        """ Writes the current executable data to a csv file.

        During the creation of the csv file, all columns (including header) are formatted based on the maximum length of all values per column.

        Parameters
        ----------
        name : str
            The ABSOLUTE path to the csv file where the data is written to.
        """
        # Skip if empty
        if self.__executable_data.empty:
            return None
        # First copy the DataFrame to a temporary object to allow formatting
        dataframe_to_write = self.__executable_data.copy()
        # Loop through all rows and call the from UserSpecifications.from_pandas to ensure that dependent data are written correctly
        for index, row in self.__executable_data.iterrows():
            user_specs = UserSpecifications.from_pandas(row)
            for spec_name, spec in user_specs.items():
                dataframe_to_write.at[index, spec_name] = spec.value
        # Apply the string operation on all entries
        dataframe_to_write = data_func.format_dataframe(dataframe_to_write, format_type=str, include_header=False)
        # Compute the maximum length of each column and format the columns based on it
        format_dict = {column: data_func.get_string_column_max_length_format(dataframe_to_write, column) for column in dataframe_to_write}
        # Format the dataframe
        dataframe_to_write = data_func.format_dataframe(dataframe_to_write, format_dict, include_header=True)
        # Write the dataframe to file
        dataframe_to_write.to_csv(name, sep=",", header=True, index=False, encoding='utf-8', na_rep="NaN")

    def add_configuration(self, xml_element: et.Element, index_based: bool = False) -> None:
        """ Adds the current executable data to a xml configuration tree.

        Parameters
        ----------
        xml_element : et.Element
            The xml tree where the data is added to. The data is added as a sub tag to this given element.
        index_based : bool, optional
            Flag whether the data is built from an index based variation or not.
        """
        # Fill the dictionary with all data. Depending on index based variations or not, all entries (index) or unique entries (non-index) are taken
        dict_to_write = {column: [] for column in self.__executable_data if column != "Executable"}
        for _, data in self.__executable_data.iterrows():
            for name, value in data.items():
                if name != "Executable":
                    if str(value) not in dict_to_write[name] or index_based:
                        dict_to_write[name].append(str(value))
        # In index based configuration get the maximum length per entry list (first of all values)
        if index_based:
            max_size = [len(value) for value in list(dict_to_write.values())[0]]
            for values in dict_to_write.values():
                max_size = [len(value) if len(value) > max_size[index] else max_size[index] for index, value in enumerate(values)]
            dict_to_write = {name: [value.ljust(max_size[index]) for index, value in enumerate(values)] for name, values in dict_to_write.items()}
        # Otherwise only adjust the values to the total maximum size
        else:
            max_size = max([len(" ".join(values)) + 1 for values in dict_to_write.values()])
            dict_to_write = {name: [(" ".join(values)).ljust(max_size)] for name, values in dict_to_write.items()}
        # Get the maximum size of all keys that are written
        max_key_size = max([len(name) + 1 for name in dict_to_write.keys()])
        # Write the values into the xml tree
        for name, value in dict_to_write.items():
            xo.modify_xml_tag(xml_element, " " * (max_key_size - len(name)) + " ".join(value), NameStyle.xml.format(name))
