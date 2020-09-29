#!/usr/bin/env python3
# Python modules
from typing import List, Tuple, Dict, Union, Optional, Any, Type, IO
import numpy as np
# alpacapy modules
from alpacapy.helper_functions.bool_type import BoolType
from alpacapy.helper_functions import check_operations as co
from alpacapy.name_style import NameStyle


class TagBase:
    """ Base class for providing an interface between InputfileTags and UserSpecificationTags.

    The TagBase class provides an interface for all specification tags (inputfile and user specification) that can be set for the Alpaca framework through
    alpacapy. It provides the functionality that both derived classes contain.
    Additional functions that differ between both are implemented in their respective class.

    Attributes
    ----------
    __value : Union[int,BoolType,str]
        The current value of the tag, by default None.
    __default_value : Union[int,BoolType,str]
        The default value used for the tag.
    __value_type : Type
        The type of the tag.
    __allowed_values : Optional[List[Union[int,BoolType,str]]]
         A list holding all allowed variables for this tag.
    """

    def __init__(self, value_type: Type, default_value: Union[int, BoolType, str],
                 allowed_values: Optional[List[Union[int, BoolType, str]]] = None) -> None:
        """Constructor

        Parameters
        ----------
        value_type : Type
            The type of the tag.
        default_value : Union[int,BoolType,str]
            The default value used for the tag.
        allowed_values : Optional[List[Union[int,BoolType,str]]], optional
            A list holding all allowed variables for this tag.
        Raises
        ------
        TypeError
            1. If float type is specified.
            2. If the type of the default value or allowed variables does not correspond to given type.
        """
        # Currently floating point types are not allowed to specify a specification (does not guarantee proper matching using == operator). Instead a string can
        # be used given the exact value, e.g. "0.000000".
        if value_type == np.dtype(float):
            raise TypeError("Currently floating values are not allowed to specify a tag. Please use string instead.")
        # Convert to wrapper bool class in case the value type is a numpy boolean
        value_type = BoolType if np.dtype(value_type) == np.dtype(bool) else np.dtype(value_type)

        # Check that the value coincides with the specified type
        if not co.check_type(default_value, value_type.type):
            raise TypeError("The provided default value does not coincide with the given type '" + str(value_type) + "'")
        # Check that all provided values coincide with the provided value type
        if allowed_values is not None and not all(co.check_type(value, value_type.type) for value in allowed_values):
            raise TypeError("Not all elements in the list of allowed variables coincide with the specified value type")

        # Assign the variables as instance variables. The type is converted into a numpy compliant data type to ensure consistency.
        self.__value = None
        self.__default_value = default_value
        self.__value_type = value_type
        self.__allowed_values = allowed_values

    def __repr__(self) -> str:
        """ Implementation of the built-in repr function """
        string = "Type           : " + str(self.__value_type) + "\n"
        string += "Value          : " + str(self.__value) + "\n"
        string += "Default value  : " + str(self.__default_value) + "\n"
        string += "Allowed values : " + str(self.__allowed_values) + "\n"
        return string

    def __str__(self) -> str:
        """ Implementation of the built-in str function """
        return type(self).__name__

    def __contains__(self, value: Union[int, float, BoolType, str]) -> bool:
        """ Implementation of the built-in contains function.

        Checks whether the provided value is in the list of allowed variables and
        corresponds to the type of the tag.

        Parameters
        ----------
        value : Union[int,float,BoolType,str]
            The value that is checked if it is contained in the list of allowed variables.
        Returns
        -------
        bool
            True if contained, False otherwise.
        """
        if self.__allowed_values is None:
            return co.check_type(value, self.__value_type.type)
        elif self.__value_type == np.dtype(BoolType):
            return value in [True, False]
        else:
            return value in self.__allowed_values

    @property
    def value(self) -> Union[int, BoolType, str]:
        """ Allows accessing the current tag value via .value. Prohibits manipulation of the variable. Therefore, the setter function must be called.
        Returns
        -------
        Union[int,BoolType,str]
            The current value of the tag.
        """
        return self.__value

    @value.setter
    def value(self, new_value: Union[int, BoolType, str]) -> None:
        """ Sets the current value to a new value.

        Parameters
        ----------
        new_value : Union[int,BoolType,str]
            the new value to be set.
        Raises
        ------
        TypeError
            If the new value does not correspond to the type.
        ValueError
            If the new value is not in the list of allowed variables.
        """
        if new_value is None:
            self.__value = None
        else:
            if not co.check_type(new_value, self.__value_type.type):
                raise TypeError("The provided value does not coincide with the type '" + str(self.__value_type) + "' of specification.")
            if self.__allowed_values is not None and new_value not in self.__allowed_values:
                raise ValueError("The given value is not allowed fof specification. Allowed variables are:\n" + str(self.__allowed_values))
            self.__value = self.__value_type.type(new_value)

    @property
    def default(self) -> Union[int, BoolType, str]:
        """ Allows accessing the default value of the tag via .value. Prohibits manipulation of the variable.
        Returns
        -------
        Union[int,BoolType,str]
            The default value of the tag.
        """
        return self.__default_value

    @property
    def allowed_values(self) -> List[Union[int, BoolType, str]]:
        """ Allows accessing the allowed variables of the tag via .allowed_variables. Prohibits manipulation of the variable.
        Returns
        -------
        List[Union[int,BoolType,str]]
            The allowed variables of the tag.
        """
        if self.__value_type == np.dtype(BoolType):
            return [True, False]
        else:
            return self.__allowed_values

    def type(self, value: Any) -> Union[str, int, BoolType]:
        """ Converts a value into the type of the tag.

        Parameters
        ----------
        value : Any
            The value that should be converted.
        Returns
        -------
        Union[str,int,BoolType]
            The converted value
        Raises
        ------
        ValueError
            If the value cannot be converted into the type of the tag.
        """
        if not co.check_type(value, self.__value_type.type):
            raise ValueError("The provided value does not coincide with the type '" + str(self.__value_type) + "' of specification.")
        return self.__value_type.type(value)

    def is_bool(self) -> bool:
        """ Checks whether the tag is a boolean tag.
        Returns
        -------
        bool
            True if it is boolean, False otherwise.
        """
        return self.__value_type == BoolType

    def is_set(self) -> bool:
        """ Checks whether the tag is currently set.
        Returns
        -------
        bool
            True if set, False otherwise.
        """
        return self.__value is not None

    def add_tag(self, name: str, use_abbreviation=True) -> str:
        """ Generates the tag with the name and value as a string.

        Adds a tag to a provided name consisting of the current value and an additional abbreviation for the name. The value is only added if present. Boolean
        values are added depending on their default values. If the boolean value corresponds to its default value no tag is added.
        The tag adding procedure can be used to name for example executables or inputfiles properly depending on their current status.
        E.g., If the maximum level of an existing inputfile is changed the add_tag function would print _ML_2. If the default value is given nothing is printed.

        Parameters
        ----------
        name : str
            The name where the current tag is added to.

        use_abbreviation : bool, optional
            Flag whether the abbreviation of the tag should be included, by default True.
        Returns
        -------
        str
            The name with the additional tag.
        """
        # In case the value of this specification is not set, do not provide any tag
        if self.__value is None:
            return ""
        if self.__value_type == np.dtype(BoolType):
            if self.__value == self.__default_value:
                return ""
            else:
                return "_No_" + name if not self.__value else "_" + name
        else:
            return "_" + NameStyle.abbreviation.format(name) + "_" + str(self.__value) if use_abbreviation else "_" + str(self.__value)
