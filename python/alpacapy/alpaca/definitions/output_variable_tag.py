#!/usr/bin/env python3
# Python modules
from typing import List, Tuple, Dict, Union, Optional, Any, Type, IO
import re
# alpacapy modules
from alpacapy.helper_functions.bool_type import BoolType
from alpacapy.helper_functions import check_operations as co
from alpacapy.name_style import NameStyle


class OutputVariableTag:
    """ Defines a single tag that can be modified as an Alpaca output variable.

    The OutputVariableTag class provides an interface for all output variables that can be set for the Alpaca framework.
    It provides the functionality to obtain tags, names in different styles and to modify the variable in the appropriate file or read its value from it.

    Attributes
    ----------
    __name : str
        The name of the variable in the specification file
    __values : List[BoolType]
        The current values of the tag (Standard, Interface, Debug), by default None.
    __default_values : List[BoolType]
        The default values used for the tag (Standard, Interface, Debug).
    """

    def __init__(self, variable_name: str, default_values: List[BoolType]) -> None:
        """Constructor

        Parameters
        ----------
        variable_name : str
            The name of the output variable to identify it uniquely in the file content.
        default_values : List[BoolType]
            The default values used for the tag (Standard, Interface, Debug).
        Raises
        ------
        TypeError
            If types of default values are not bools.
        ValueError
            1. If the variable name is not a string or is empty
            2. If the number of default values is not three.
        """
        # Assign value type
        value_type = BoolType

        # Check that the value coincides with the specified type and are three elements
        if len(default_values) != 3:
            raise ValueError("The number of default values for output variables must be three.")
        for default_value in default_values:
            if not co.check_type(default_value, value_type.type):
                raise TypeError("The provided default value does not coincide with the given type '" + str(value_type) + "'")

        # Check that the variable name is not empty
        if not isinstance(variable_name, str) or not variable_name:
            raise ValueError("The output variable name must be of type string and non-empty")

        # Assign the variables as instance variables. The type is converted into a numpy compliant data type to ensure consistency.
        self.__name = variable_name
        self.__value_type = BoolType
        self.__values = None
        self.__default_values = [self.__value_type.type(value) for value in default_values]

    def __repr__(self) -> str:
        """ Implementation of the built-in repr function """
        string = "Type               : " + str(self.__value_type) + "\n"
        string += "Values             : " + str(self.__values) + "\n"
        string += "Default values     : " + str(self.__default_values) + "\n"
        return string

    def __str__(self) -> str:
        """ Implementation of the built-in str function """
        return type(self).__name__

    @property
    def values(self) -> Union[int, BoolType, str]:
        """ Allows accessing the current tag value via .value. Prohibits manipulation of the variable. Therefore, the setter function must be called.
        Returns
        -------
        Union[int,BoolType,str]
            The current value of the tag.
        """
        return self.__values

    @values.setter
    def values(self, new_values: List[BoolType]) -> None:
        """ Sets the current value to a new value.

        Parameters
        ----------
        new_value : Union[int,BoolType,str]
            the new value to be set.
        Raises
        ------
        TypeError
            If the new value does not correspond to the type.
        """
        if new_values is None:
            self.__values = None
        else:
            if len(new_values) != 3:
                raise ValueError("The number of values for output variables must be three.")
            for value in new_values:
                if not co.check_type(value, self.__value_type.type):
                    raise TypeError("The provided value does not coincide with the type '" + str(self.__value_type) + "' of output variable.")
            self.__values = [self.__value_type.type(value) for value in new_values]

    @property
    def default(self) -> Tuple[BoolType]:
        """ Allows accessing the default values of the tag via .default. Prohibits manipulation of the variable.
        Returns
        -------
        List[BoolType]
            The default values of the tag for standard, interface and debug output.
        """
        return self.__default_values

    def modify_file(self, file_content: str, use_default: bool = False) -> str:
        """ Modifies the output variable in the provided file content.

        Parameters
        ----------
        file_content : str
            The content of the specification file. The file content must be fully loaded into cache before calling this function.
        use_default : bool, optional
            Flag whether the specification default value should be used or not, by default False.
        Returns
        -------
        str
            The modified file_content.
        """
        if self.values is None and not use_default:
            return file_content
        # Create the c++ compliant true or false
        values = self.__default_values if use_default else self.__values
        values = ["true" if value else "false" for value in values]
        return re.sub("( " + self.__name + " *)=( *{)(.+?),(.+?),(.+?)};",
                      r"\g<1>=\g<2> " + values[0] + ", " + values[1] + ", " + values[2] + " };", file_content)

    def read_from_file(self, file_content: str) -> None:
        """ Reads the data for this specification from the specification file.

        Parameters
        ----------
        file_content : str
            The content of the specification file. The file content must be fully loaded into cache before calling this function.
        """
        # Search the appropriate tags for this variable
        read_value = re.search("( " + self.__name + " *)=( *{)(.+?),(.+?),(.+?)};", file_content)
        if read_value is None:
            ValueError("Error reading the output variable. Maybe the wrong variable name is specificed.")
        self.values = [read_value.group(index).strip() for index in [3, 4, 5]]

    def add_tag(self, name: str, use_abbreviation=True) -> str:
        """ Generates the tag with the name and value as a string.

        This functions is required to align the OutputVariable to the other specifications.

        Parameters
        ----------
        name : str
            The name where the current tag is added to.

        use_abbreviation : bool, optional
            Flag whether the abbreviation of the tag should be included, by default True.
        Returns
        -------
        str
            Always an empty string since output
        """
        if not any(self.__values):
            return ""
        else:
            tag = "_" + NameStyle.abbreviation.format(name) if use_abbreviation else "_" + name
            tag += "_"
            for type, value in zip(["S", "I", "D"], self.__values):
                if value:
                    tag += type
            return tag
