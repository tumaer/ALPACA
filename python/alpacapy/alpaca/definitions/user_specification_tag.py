#!/usr/bin/env python3
# Python modules
from typing import List, Tuple, Dict, Union, Optional, Any, Type, IO
import numpy as np
import re
# alpacapy modules
from alpacapy.alpaca.definitions.tag_base import TagBase
from alpacapy.alpaca.definitions.user_specification_file import UserSpecificationFile
from alpacapy.helper_functions.bool_type import BoolType


class UserSpecificationTag(TagBase):
    """ Defines a single tag that can be modified as an Alpaca user specifications.

    The UserSpecificationTag class provides an interface for all user specification variables that can be set for the Alpaca framework.
    It provides the functionality to obtain tags, names in different styles and to modify the variable in the appropriate file or read its value from it.

    Parameters
    ----------
    TagBase :
        Base class to provide generic interface for most functions to InputfileTag.
    Attributes
    ----------
    __specification_filename : str
        The filename where the user specification can be found.
    __specification_file_tag : str
        The tag that clearly identifies the specification in the file (often with a trailing underscore).
    __specification_file_prefix : str
        Additional prefix that must be used to replace the value in the file, by default None (e.g., a Namespace::).
    """

    def __init__(self, value_type: Type, default_value: Union[int, BoolType, str],
                 specification_file: UserSpecificationFile, specification_file_tag: str, specification_file_prefix: Optional[str] = None,
                 allowed_values: Optional[List[Union[int, BoolType, str]]] = None, no_tag: Optional[bool] = False) -> None:
        """Constructor

        Parameters
        ----------
        value_type : Type
            The type of the value that corresponds to the specification (bool,int,float,str).
        default_value : Union[int,BoolType,str]
            The default value used for the specification.
        specification_file : UserSpecificationFile
            The file in the Alpaca/src/user_specification folder where this specification is defined.
        specification_file_tag : str
            The tag that clearly identifies the specification in the file (often with a trailing underscore).
        specification_file_prefix : Optional[str], optional
            Additional prefix that must be used to replace the value in the file, by default None (e.g., a Namespace::).
        allowed_values : Optional[List[Union[int,float,BoolType,str]]], optional
            A list holding all allowed values for this specification, by default None
        no_tag : Optional[bool]
            Flag whether tags for this specification should be printed or not, by default False.
        Raises
        ------
        ValueError
            If the provided user specification file is not in the list of allowed file names.
        """
        # Consistency check of variables
        if specification_file not in UserSpecificationFile:
            raise ValueError("The provided user specification file is not contained in the allowed variables of UserSpecificationFile")

        # Call the base class constructor with appropriate values
        super().__init__(value_type, default_value, allowed_values, no_tag)
        # Assign all class instance variables
        self.__specification_filename = specification_file
        self.__specification_file_tag = specification_file_tag
        self.__specification_file_prefix = specification_file_prefix

    def __repr__(self) -> str:
        """ Implementation of the built-in repr function """
        string = super().__repr__()
        string += "Specification file       : " + str(self.__specification_filename) + "\n"
        string += "Specification file tag   : " + str(self.__specification_file_tag) + "\n"
        string += "Specification file prefix: " + str(self.__specification_file_prefix) + "\n"
        return string

    @property
    def filename(self) -> str:
        """ Allows accessing the specification filename via .filename.

        Prohibits manipulation of the variable.
        Returns
        -------
        str
            The filename attribute of the class.
        """
        return self.__specification_filename

    def modify_file(self, file_content: str, use_default: bool = False) -> str:
        """ Modifies the specification in the provided file content.

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
        if self.value is None and not use_default:
            return file_content
        value_to_set = self.default if use_default else self.value
        if self.is_bool():
            value_to_set = "true" if self.value else "false"
        prefix = "" if self.__specification_file_prefix is None else self.__specification_file_prefix
        # The substitution is important to maintain the leading white space before and after the = operator. Otherwise it cannot guaranteed that only the
        # tag is modified.
        return re.sub(" (" + self.__specification_file_tag + " *)= (.+?);",
                      r" \g<1>= " + prefix + str(value_to_set) + ";", file_content)

    def read_from_file(self, file_content: str) -> None:
        """ Reads the data for this specification from the specification file.

        Parameters
        ----------
        file_content : str
            The content of the specification file. The file content must be fully loaded into cache before calling this function.
        """
        # Specify the different groups (tag everythin in front of the = sign, prefix and suffix to specify everything before and after the value to read)
        tag_group = "(.+?" + self.__specification_file_tag + ".+?)"
        prefix_group = "" if self.__specification_file_prefix is None else "(.+?" + self.__specification_file_prefix + ")"
        # Set the internal value via the set function (makes consistency checks)
        read_value = re.search(tag_group + "=" + prefix_group + "(.+?);", file_content)
        if read_value is None:
            ValueError("Error reading the specification variable. Maybe the wrong tag is specificed.")
        self.value = read_value.group(2 if not prefix_group else 3).strip()
